'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers coverage metrics data and prepares it for the output report.
It relies on the following official sources:
https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/CoverageMetricsReport_fDG.htm
https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/QCMetricsCoverageReports_fDG_dtSW.htm
Section: DRAGEN DNA Pipeline, Subsection: QC Metrics and Coverage/Callability Reports
https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf

The following sentence from the User-Guide deserves attention:
Additional coverage metrics can be enabled, and additional coverage regions can be specified.
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

import logging
import re
from collections import defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

from .utils import check_duplicate_samples, clean_headers, make_log_report, order_headers

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "Coverage"


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The following code section provides a clean and simple tool for setting configurations:
https://github.com/MultiQC/MultiQC/blob/main/docs/plots.md#creating-a-table

The SINGLE_HEADER is used to define the common settings for all coverage headers.
'''
SINGLE_HEADER = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "title": None,  # Short title, table column title
    "description": None,  # Longer description, goes in mouse hover text
    "max": None,  # Minimum value in range, for bar / colour coding
    "min": 0,  # Maximum value in range, for bar / colour coding
    "ceiling": None,  # Maximum value for automatic bar limit
    "floor": None,  # Minimum value for automatic bar limit
    "minrange": None,  # Minimum range for automatic bar
    "scale": "GnBu",  # Colour scale for colour coding. False to disable.
    "bgcols": None,  # Dict with values: background colours for categorical data.
    "colour": "15, 150, 255",  # Colour for column grouping
    "suffix": None,  # Suffix for value (eg. "%")
    "format": "{:,.2f}",  # Value format string - MultiQC default 1 decimal places
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values, special case, see below
    "hidden": True,  # Set to True to hide the column in the general table on page load.
    "hidden_own": True,  # For non-general plots in own coverage sections.
    "exclude": True,  # Exclude all headers from the general html table.
}
'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER. These are just extensions,
but can be used as if they were a part of the SINGLE_HEADER. Extra configurations are not
set by default. Only the SINGLE_HEADER is used as a basis for a header. But you can always
append an extra config(eg hidden_own, exclude) to the SINGLE_HEADER to overcome this issue.

Hints and rules:
- "hidden" is used for own coverage section if "hidden_own" is not provided.

- "exclude" and "exclude_own" may be useful in cases where html size matters.
  Eg. many samples, several regions but only few metrics of interest.

- "order_priority" shall be used only in the METRICS.
  Details are in the corresponding comment below.
'''
EXTRA_HEADER = {
    "hidden_own": None,  # For non-general plots in the coverage section.
    "exclude": None,  # True to exclude metric from the general html table.
    "exclude_own": None,  # True to exclude from own coverage section html table.
    "order_priority": None,  # Used to specify columns' order in all tables.
}
'''"""""""""""""""""""""""""""""""""""""""""""""""""""""
The composition of each region string is very sensitive.
It has to coincide with the region format in input-CSVs.
Letters must be lower case. No padding(extra spaces) on
the left/right sides. Words must be separated from each
other with a single-space string, that is " ".
'''
WGS = "genome"
QC = "qc coverage region"
BED = "target region"
# The desired regions can be added:
# REGION = "some region"

"""
The REGIONS holds the variables. You can add common region's settings here.
These configs have higher priority than those in the SINGLE_HEADER.
They also overwrite the automatically created ones.
REGIONS is also used to detect actual regions in real data and allow
the special metrics "aligned bases/reads" to have region-specific settings.
So it is a good idea to include all known regions in the dictionary.
"""
REGIONS = {
    WGS: {
        # "hidden": False,
        # "scale": "Blues",
    },
    QC: {
        # "exclude_own": True,
        # "scale": "Greens",
    },
    BED: {
        # "scale": "Reds",
    },
}

'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The METRICS is the container for abstract representations of the standard metrics.
They are designed to be seen as simple identifiers of the well known metrics.
You can set any possible real/virtual configuration defined in the SINGLE_HEADER.

Some rules:
- You can use the single regions defined above as keys to set
  region-specific configurations for a certain metric.

- "modify" must conform to the following general structure:
  lambda x: (arbitrary expression) if isinstance(x, str) else (math expression)
  The x is either number or string. The latter is some value, which could not be
  converted to int or float. NA for instance. Such values can be modified in the
  (arbitrary expression) if needed. If not then just return the input.

- "order_priority" can be used to define an insertion order of headers.
  Only int and float numbers are valid. The config is used for own and
  general tables.

- The descending precedence of configurations:
  1. Region-specific metric
  2. General metric
  3. General region
  4. Automatically defined settings
  5. Basic configs in the SINGLE_HEADER
""""""""""""""""""""""""""""""""""""'''

""" Extra Info:

- If you want to add new metrics, then you could try the following:
  * All characters must be in lower case.
  * No extra padded space on the left and right sides.
  * If tokens (eg words, operators) are separated from each other by a sequence
    of spaces, replace each such sequence with just one single-space string " ".
  * If it still does not work, then you might need to take a look
    at the "make_consistent_metric" and maybe at the "make_user_configs".

- If you want to use "greater than or equal to" sign, then take a look at these:
    ">="
    &ge;
    ⩾
    ≥
    &#8925;  or  &#x22DD;
    &#8927;  or  &#x22DF;
    &#8829;  or  &#x227D;

  New characters can be added as well.
  Some of them may be more pleasing to your eye than others.

- There is a little harmless bug in the process of MultiQC table creation.
  It might be of interest if you are using the "exclude" config. So,
  if there is only one column in the table, and the "hidden" set to True,
  then the column !AND! the "Configure Columns" button won't be shown.
  Unfortunately there is no simple user-friendly way to extract column's data.
  The column is actually in the html table, but the "display" attribute is set
  to none. "Configure Columns" button is not included in table's div block.

- The "cond_formatting_colours" along with "cond_formatting_rules" can be
  used to put extra colour on values in html tables. The former must be a
  list of dictionaries, where each dictionary shall be a single definition
  of some colour. The latter is a dictionary with colours. Each colour's
  value is a list, which holds dictionaries with ecomparison rules, which
  can be found in the little code block below. Please notice that an insertion
  order of colours in the "cond_formatting_colours" plays an enormous role,
  because the last chosen colour will be used in html.

  The code can be modified in the future, but probably not dramatically.
  So the following description could be useful. At least for some time.

  # This python code is a simplified part of the multiqc/plots/table.py
  # It shows the logical process of defining a colour on some value.
  VAL # The actual single input value. Most likely int or float.
      # It's been modified by the "modify" lambda if it was defined.
  APPLY # Applies colour on the VAL. Abstract, non-existent function.

  for colour in cond_formatting_colours:      # User + MultiQC colours.
    match = False
    for cmp in cond_formatting_rules[colour]: # User + MultiQC rules.
      try:
        # Each comparison should be a dict with single key: val
        if "s_eq" in cmp and str(cmp["s_eq"]).lower() == str(VAL).lower():
          match = True
        if "s_contains" in cmp and str(cmp["s_contains"]).lower() in str(VAL).lower():
          match = True
        if "s_ne" in cmp and str(cmp["s_ne"]).lower() != str(VAL).lower():
          match = True
        if "eq" in cmp and float(cmp["eq"]) == float(VAL): match = True
        if "ne" in cmp and float(cmp["ne"]) != float(VAL): match = True
        if "gt" in cmp and float(cmp["gt"]) <  float(VAL): match = True
        if "lt" in cmp and float(cmp["lt"]) >  float(VAL): match = True
      except:
        logger.warning("Not able to apply table conditional formatting")
    if match:
        APPLY(colour)

- Multiline strings are not used in "description" to guarantee that indents
  won't be included. Though browsers reduce such chains to just one space,
  it still feels safer to avoid them by concatenating single strings.


Ideas:

- "format" could be automatically updated depending on the input value.
  if int -> "{:,.0f}
  elif float -> "{:,.Nf}", where N is the number of after-point digits.
  else -> "{:,.1f}, which is default by the MultiQC.

Final notes:

- Most configs below are for demonstration purposes only.
  You can adapt them according to your personal taste/needs.

- You have the final word in defining settings.
  The code trusts you and does not check the chosen values.
"""

# These might be useful for checking results.
AUTO_CONFIGS_ENABLED = True
USER_CONFIGS_ENABLED = True

V2 = " pct"  # Used to create a unique ID if a metric has two values.

METRICS = {
    "aligned reads": {
        "order_priority": 0,
        "hidden_own": False,
        "colour": "255, 0, 0",
    },
    "aligned reads in region": {
        "order_priority": 0.1,
        "scale": "RdGy",
        "colour": "255, 0, 0",
    },
    "aligned reads in region" + V2: {
        "order_priority": 0.2,
        "max": 100,
        "suffix": "%",
        "scale": "Purples",
        "colour": "255, 0, 0",
        "exclude": False,
    },
    "aligned bases": {
        "order_priority": 0.3,
        "hidden_own": False,
        "scale": "Greens",
        "colour": "0, 0, 255",
    },
    "aligned bases in region": {
        "order_priority": 0.4,
        "scale": "Reds",
        "colour": "0, 0, 255",
        WGS: {
            "scale": "Purples",
        },
        QC: {
            "scale": "Oranges",
        },
    },
    "aligned bases in region" + V2: {
        "order_priority": 0.5,
        "max": 100,
        "colour": "0, 0, 255",
        "scale": "Greens",
        "suffix": "%",
        "bgcols": {"NA": "#00FFFF"},
        "exclude": False,
    },
    "average chr x coverage over region": {
        "order_priority": 1.0,
        "title": "X cov",
        "suffix": " x",
        "colour": "179,179,50",  # Olive
    },
    "average chr y coverage over region": {
        "order_priority": 1.1,
        "title": "Y cov",
        "suffix": " x",
        "scale": "YlOrRd",
        "colour": "179,179,50",  # Olive
    },
    "xavgcov/yavgcov ratio over region": {
        "order_priority": 1.2,
        "title": "XAvgCov/YAvgCov",
    },
    "xavgcov/autosomalavgcov ratio over region": {
        "order_priority": 1.3,
        "title": "XAvgCov/AutAvgCov",
    },
    "yavgcov/autosomalavgcov ratio over region": {
        "order_priority": 1.4,
        "title": "YAvgCov/AutAvgCov",
    },
    "average alignment coverage over region": {
        "order_priority": 3,
        "exclude": False,
        "hidden_own": False,
        "hidden": False,
        "title": "Depth",
        "scale": "BrBG",
        "colour": "55,126,184",  # Blue
    },
    "average mitochondrial coverage over region": {
        "order_priority": 4,
        "title": "MT cov",
        "scale": "Purples",
        "suffix": " x",
        "colour": "152,78,163",  # Purple
    },
    "average autosomal coverage over region": {
        "order_priority": 6,
        "title": "Mean aut cov",
        "suffix": " x",
        "scale": "Reds",
        "bgcols": {"NA": "#FF00FF"},
        WGS: {
            "description": "Average autosomal coverage over genome. Calculated"
            " as the number of bases that aligned to the autosomal loci in genome"
            " divided by the total number of loci in the autosomal loci in"
            " genome. If there are no autosomes in the reference genome, or"
            " the region does not intersect autosomes, this metric shows as NA.",
        },
    },
    "median autosomal coverage over region": {
        "order_priority": 7,
        "title": "Med aut cov",
        "suffix": " x",
        "scale": "Greens",
        "exclude": False,
    },
    "mean/median autosomal coverage ratio over region": {
        "order_priority": 7.1,
        "title": "Mean/med autosomal coverage",
        "hidden_own": False,
        "suffix": " x",
    },
    # This metric can handle any float.
    # Title and description shall not be set here, they are automated.
    # But if you want to handle a certain number manually, then you can
    # add a float-string key with desired configs in the "extra".
    "uniformity of coverage (pct > d.d*mean) over region": {
        "suffix": "%",
        "colour": "77,175,74",  # Green
        "exclude": False,
        "extra": {
            # "Uniformity of coverage (PCT > 0.2*mean) over region" metric.
            "0.2": {
                "scale": "Oranges",
                "order_priority": 5.2,
                "hidden_own": False,
                # "Uniformity of coverage (PCT > 0.2*mean) over QC coverage region" metric.
                QC: {
                    "description": "Percentage of sites with coverage greater than"
                    " 20% of the mean coverage in QC coverage region. Demonstrates"
                    " the uniformity of coverage, the higher the better.",
                },
                # "Uniformity of coverage (PCT > 0.2*mean) over genome" metric.
                WGS: {
                    "description": "Percentage of sites with coverage greater"
                    " than 20% of the mean coverage in genome. Demonstrates"
                    " the uniformity of coverage, the higher the better.",
                },
            },
            # "Uniformity of coverage (PCT > 0.4*mean) over region" metric.
            "0.4": {
                "order_priority": 5.4,
                QC: {
                    "description": "Percentage of sites with coverage greater"
                    " than 40% of the mean coverage in QC coverage region.",
                },
            },
        },
    },
    # The following two metrics are the most tedious. The creation of the title
    # and description are automated. If you want to target some specific case,
    # then you can specify a certain combination of (ix, inf) or (ix, jx).
    # Just write a group in the "extra", you can consider it as a separate
    # metric, so you can also add general or region-specific configs there.
    "pct of region with coverage [ix:inf)": {
        "max": 100,
        "scale": "Purples",
        "suffix": "%",
        "colour": "228,26,28",  # Red
        "exclude": False,
        WGS: {
            # "scale": "Reds",
        },
        QC: {
            # "scale": "Oranges",
        },
        "extra": {
            # "PCT of region with coverage [0x: inf)" metric.
            ("0x", "inf"): {
                "title": "≥0x",
                # Simple example of colors with rules.
                # "cond_formatting_rules": {
                #    "red_bad": [{"s_contains": ""}], # Each value is red by default.
                #    "green_good": [{"eq": 100.0}],   # If equal to 100, then green.
                # },
                # "cond_formatting_colours": [
                #    {"red_bad": "#FF0000"},
                #    {"green_good": "#00FF00"},
                # ],
                # "PCT of genome with coverage [0x: inf)" metric.
                WGS: {
                    "description": "Percentage of sites in genome with any coverage.",
                },
                # "PCT of QC coverage region with coverage [0x: inf)" metric.
                QC: {
                    "description": "Percentage of sites in QC coverage region with any coverage.",
                },
                # "PCT of target region with coverage [0x: inf)" metric.
                BED: {
                    "description": "Percentage of sites in target region with any coverage.",
                },
            },
            # "PCT of region with coverage [1x: inf)" metric.
            ("1x", "inf"): {
                # "title": "Some title",
                # "description": "Some description",
            },
            # "PCT of region with coverage [100x: inf)" metric.
            ("100x", "inf"): {
                # "title": "Another title",
                # "description": "Some description",
            },
        },
    },
    "pct of region with coverage [ix:jx)": {
        "max": 100,
        "scale": "Reds",
        "suffix": "%",
        # "exclude_own": True,
        "extra": {
            ("0x", "1x"): {
                "title": "0x",
                WGS: {
                    "description": "Percentage of sites in genome with no coverage.",
                },
                QC: {
                    "description": "Percentage of sites in QC coverage region with no coverage.",
                },
                BED: {
                    "description": "Percentage of sites in target region with no coverage.",
                },
            },
        },
    },
}

# The TABLE_CONFIG defines common configs for all whole tables.
TABLE_CONFIG = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "id": None,  # ID used for the table
    "table_title": None,  # Title of the table. Used in the column config modal
    "save_file": False,  # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sort_rows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_violin": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}
# Below are region-specific configs with higher priority.
REGION_TABLE_CONFIG = {
    WGS: {
        "table_title": "WGS",
        # "no_violin": True,
        # "save_file": True,
        # "raw_data_fn": "Genome_raw_file",
    },
    QC: {
        "table_title": "QC Coverage Region",
    },
    BED: {
        "table_title": "Target Bed",
    },
}


##**********************************************************************************************##

##                                  The actual code starts here

##**********************************************************************************************##


class DragenCoverageMetrics(BaseMultiqcModule):
    """Public members of the DragenCoverageMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined.
    Other methods can be added as well to provide extra features (eg module interface, JSON).
    """

    def add_coverage_metrics(self, overall_mean_cov_data):
        """The main function of the dragen coverage metrics module.
        The public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns a set with sample names."""

        cov_data = defaultdict(dict)

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        # Used to properly extract data from *_overall_mean_cov.csv files.
        match_overall_mean_cov = defaultdict(dict)

        # Metric IDs with original consistent metric strings.
        all_metrics = {}

        for file in self.find_log_files("dragen/coverage_metrics"):
            out = coverage_parser(file)
            if out["success"]:
                self.add_data_source(file, section="coverage_metrics")

                original_sample = out["sample_name"]
                cleaned_sample = self.clean_s_name(original_sample, file)
                phenotype = out["phenotype"]

                cov_data[cleaned_sample][phenotype] = out["data"]  # Add/overwrite the sample.

                all_samples[cleaned_sample].append(file)
                match_overall_mean_cov[cleaned_sample][phenotype] = (original_sample, file["root"])
                all_metrics.update(out["metric_IDs_with_original_names"])

                # Superfluous function call to confirm that it is used in this module
                # Replace None with actual version if it is available
                self.add_software_version(None, cleaned_sample)

        cov_data = self.ignore_samples(cov_data)
        if not cov_data:
            return set()

        cov_headers = make_coverage_headers(all_metrics)

        # Check samples for duplicates.
        check_duplicate_samples(all_samples, log)

        # Report found info/warnings/errors, which were collected while
        # calling the coverage_parser and constructing cov_headers.
        make_log_report(log_data, log, "dragen/coverage_metrics")

        # Write data into the general table.
        gen_data, gen_headers = make_general_stats(cov_data, cov_headers)
        self.general_stats_addcols(gen_data, gen_headers, namespace=NAMESPACE)

        # Write data to file.
        out_data = make_data_for_txt_reports(cov_data, all_metrics)
        for phenotype in out_data:
            self.write_data_file(out_data[phenotype], phenotype)

        # Extract coverage bed/target bed/wgs from _overall_mean_cov.csv files.
        # And prepare <coverage-region-prefix>-specific texts.
        bed_texts = make_bed_texts(overall_mean_cov_data, match_overall_mean_cov)
        coverage_sections = make_cov_sections(cov_data, cov_headers, bed_texts)

        for cov_section in coverage_sections:
            self.add_section(
                name=cov_section["name"],
                anchor=cov_section["anchor"],
                description=cov_section["description"],
                helptext=cov_section["helptext"],
                plot=cov_section["plot"],
            )
        return cov_data.keys()


def make_data_for_txt_reports(coverage_data, metrics):
    """Prepare phenotype-specific data for text reports."""

    data = defaultdict(lambda: defaultdict(dict))
    V2_len = len(V2)  # Used to "cut off" the V2 from the metric.
    for sample in coverage_data:
        for phenotype in coverage_data[sample]:
            # Replace any sequence of spaces/hyphens/dots/underscores by single underscore.
            new_phenotype = re.sub(r"(\s|-|\.|_)+", "_", phenotype)
            # Transform phenotype as in the previous code version.
            if new_phenotype == "wgs":
                new_phenotype += "_cov_metrics"
            elif "qc_coverage_region" in new_phenotype:
                new_phenotype = new_phenotype.replace("qc_coverage_region", "qc_region") + "_coverage_metrics"
            else:
                new_phenotype += "_coverage_metrics"
            new_phenotype = "dragen_" + new_phenotype

            for metric, value in coverage_data[sample][phenotype]["data"].items():
                # If metric is for the second value, then it has V2 at the end.
                # "Cut" it, extract original consistent metric and append V2.
                original_metric = metrics[metric] if metric in metrics else metrics[metric[:-V2_len]] + V2
                data[new_phenotype][sample][original_metric] = value
    return data


def make_bed_texts(overall_mean, coverage_data):
    """Matches _overall_mean_cov.csv to the corresponding _coverage_metrics.csv
    Extracts coverage bed/target bed/wgs/file names and creates a text for each
    section (phenotype)."""

    # Each sample.phenotype can have at most 1 corresponding overall_mean_cov.csv file.
    # If it has it, then append the sample to the sources_matched. If not, then
    # append to the sources_not_matched, in order to properly deliver info to users.
    sources_matched = defaultdict(list)
    sources_not_matched = defaultdict(list)
    phenotypes = set()  # Stores all found phenotypes.

    for sample in coverage_data:
        for phenotype in coverage_data[sample]:
            phenotypes.add(phenotype)
            orig_sample, root = coverage_data[sample][phenotype]
            # Check if that file is present in the overall_mean_cov data.
            if (
                overall_mean
                and root in overall_mean
                and orig_sample in overall_mean[root]
                and phenotype in overall_mean[root][orig_sample]
            ):
                # data has 2 keys: "source_file" and "value"
                data = overall_mean[root][orig_sample][phenotype]
                sources_matched[phenotype].append((sample, data["source_file"]))
            else:
                sources_not_matched[phenotype].append(sample)

    def extract_source(source_file):
        if "/" in source_file:
            return source_file.split("/")[-1]

        elif "\\" in source_file:
            return source_file.split("\\")[-1]

        # Probably just a file without path.
        else:
            return source_file

    texts_for_sections = defaultdict(dict)
    # Now go through the collected phenotypes.
    for phenotype in phenotypes:
        # This separation is used to handle the case when input directory has many different
        # sets with samples. So the shorter text goes into description, the longer into 'Help'.
        text_help_dropdown = ""
        text_description = "Information about source files can be found in the `Help` drop-list."

        # Does at least 1 sample has a source file?
        if phenotype in sources_matched:
            # There are probably not too many sources according to examined data.
            # Restructure first, so each source references list with samples.
            bed_sources = defaultdict(list)
            for sample, source in sources_matched[phenotype]:
                bed_sources[source].append(sample)

            # Maybe not all samples were matched.
            if phenotype in sources_not_matched:
                # First tell about matched ones.
                for source in bed_sources:
                    text_help_dropdown += (
                        "\n\nThe following samples are based on the " + extract_source(source) + ":\n\n"
                    )
                    for sample in bed_sources[source]:
                        text_help_dropdown += sample + ", "
                    text_help_dropdown = text_help_dropdown[:-2] + "."  # Get rid of the last ", " and append a "."

                # Also tell about unmatched samples.
                text_help_dropdown += "\n\nThe following samples do not have associated source files:\n\n"
                for sample in sources_not_matched[phenotype]:
                    text_help_dropdown += sample + ", "
                text_help_dropdown = text_help_dropdown[:-2] + "."  # Get rid of the last ", " and append a "."

            # Each sample has some corresponding overall_mean_cov file.
            else:
                # Just 1 source file for all samples?
                if len(bed_sources) == 1:
                    text_description = "All samples are based on the " + extract_source(list(bed_sources)[0]) + "."

                # There are at least 2 source files.
                else:
                    for source in bed_sources:
                        text_help_dropdown += (
                            "\n\nThe following samples are based on the " + extract_source(source) + ":\n\n"
                        )
                        for sample in bed_sources[source]:
                            text_help_dropdown += sample + ", "
                        text_help_dropdown = text_help_dropdown[:-2] + "."  # Get rid of the last ", " and append "."

        # Otherwise not a single sample has a source from overall_mean_cov.csv.
        else:
            # Try to make the text more specific.
            if re.search("^wgs$", phenotype, re.IGNORECASE):
                text_description = "No wgs source file found."
            else:
                text_description = "No 'coverage bed/target bed' source files found."

        texts_for_sections[phenotype]["description"] = text_description
        texts_for_sections[phenotype]["helptext"] = text_help_dropdown

    return texts_for_sections


def create_table_handlers():
    """Codeblock for handling tables.
    Output 2 closures:
    * make_general_stats
    * make_own_coverage_sections"""

    # Regexes for the well-known standard phenotypes.
    RGX_WGS_PHENOTYPE = re.compile(r"^wgs$", re.IGNORECASE)
    RGX_QC_PHENOTYPE = re.compile(r"^qc-coverage-region-(?P<number>\d+)$", re.IGNORECASE)
    RGX_BED_PHENOTYPE = re.compile(r"^target_bed$", re.IGNORECASE)

    def improve_gen_phenotype(phenotype):
        """Improve phenotype's appearance for columns' titles in the general table."""

        pheno_match = RGX_QC_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "QC-" + pheno_match["number"]

        pheno_match = RGX_WGS_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "WGS"

        pheno_match = RGX_BED_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "BED"

        return phenotype

    def make_general_stats(coverage_data, coverage_headers):
        """Prepare data and headers for the general table."""

        gen_data = defaultdict(dict)
        gen_headers = {}
        for sample in coverage_data:
            for phenotype in coverage_data[sample]:
                data = coverage_data[sample][phenotype]["data"]
                region = coverage_data[sample][phenotype]["region"]
                for metric in data:
                    _metric = metric
                    if metric == "aligned reads" or metric == "aligned bases":
                        if metric + region in coverage_headers:
                            _metric += region
                    # Data and the corresponding header are included in the report,
                    # only if "exclude" is not present or False/False-equivalent.
                    if not ("exclude" in coverage_headers[_metric] and coverage_headers[_metric]["exclude"]):
                        # Make exclusive metric ID.
                        # Please notice that special signs (eg "]") are
                        # excluded when HTML IDs are created. So, for example:
                        # PCT of region with coverage [10x: 50x)
                        # PCT of region with coverage [10x: 50x]
                        # will both reference the same HTML ID.
                        # Irrelevant in this module, but may not in general.
                        m_id = re.sub(r"(\s|-|\.|_)+", " ", phenotype + "_" + metric)
                        gen_data[sample][m_id] = data[metric]
                        gen_headers[m_id] = coverage_headers[_metric].copy()
                        del gen_headers[m_id]["colour"]
                        """
                        Some modifications are necessary to improve informativeness
                        of the general table, because several/many familiar tables
                        can be combined and inserted into it. HTML entities are also
                        deleted, because the title shall become wide enough after
                        concatenating the phenotype.
                        Check if "title" is present for safety, because
                        it can be set to None in the METRICS. Silly, i know.
                        """
                        if "title" in gen_headers[m_id]:
                            gen_headers[m_id]["title"] = (
                                gen_headers[m_id]["title"] + " " + improve_gen_phenotype(phenotype)
                            )

        return gen_data, clean_headers(order_headers(gen_headers))

    def make_section_name(phenotype):
        """Improve phenotype's appearance for non-general sections' titles."""

        pheno_match = RGX_QC_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "QC Region " + pheno_match["number"] + " Coverage Metrics"

        pheno_match = RGX_WGS_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "WGS Coverage Metrics"

        pheno_match = RGX_BED_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "Target Bed Coverage Metrics"

        # Try to make it better looking a little bit.
        phenotype = re.sub(r"(\s|-|\.|_)+", " ", phenotype).strip()
        phenotype = phenotype[0].capitalize() + phenotype[1:]
        return phenotype + " Coverage Metrics"

    # Regexes for the well-known regions:
    RGX_WGS_REGION = re.compile("^genome$")
    RGX_QC_REGION = re.compile("^qc coverage region$")
    RGX_BED_REGION = re.compile("^target region$")

    def improve_region(region):
        """Modify the extracted lower-case region."""

        if RGX_QC_REGION.search(region):
            return "QC coverage region"
        if RGX_WGS_REGION.search(region):
            return "genome"
        if RGX_BED_REGION.search(region):
            return "target Bed"

        return region

    def make_own_coverage_sections(coverage_data, coverage_headers, bed_texts):
        """Create non-general phenotype-specific sections."""

        plots = defaultdict(lambda: defaultdict(dict))
        # There is no guarantee, that all files would have the same region.
        regions = defaultdict(set)
        for sample in coverage_data:
            for phenotype in coverage_data[sample]:
                real_data = coverage_data[sample][phenotype]["data"]
                region = coverage_data[sample][phenotype]["region"]
                data = {sample: {}}
                headers = {}
                for metric in real_data:
                    _metric = metric
                    if metric == "aligned reads" or metric == "aligned bases":
                        if metric + region in coverage_headers:
                            _metric += region
                    if not ("exclude_own" in coverage_headers[_metric] and coverage_headers[_metric]["exclude_own"]):
                        m_id = re.sub(r"(\s|-|\.|_)+", " ", phenotype + "_" + metric)
                        data[sample][m_id] = real_data[metric]
                        headers[m_id] = coverage_headers[_metric].copy()

                        if "hidden_own" in headers[m_id]:
                            headers[m_id]["hidden"] = headers[m_id]["hidden_own"]

                plots[phenotype]["data"].update(data)
                plots[phenotype]["headers"].update(headers)
                regions[phenotype].add(region)

                plots[phenotype]["config"] = {config: val for config, val in TABLE_CONFIG.items() if val is not None}
                if region in REGION_TABLE_CONFIG and REGION_TABLE_CONFIG[region]:
                    plots[phenotype]["config"].update(REGION_TABLE_CONFIG[region])

        sections = []
        helptext = (
            "The following criteria are used when calculating coverage:\n\n"
            "* Duplicate reads and clipped bases are ignored.\n\n"
            "* DRAGEN V3.4 - 3.7: Only reads with `MAPQ` > `min MAPQ` and bases with `BQ` > `min BQ` are considered\n\n"
            "* DRAGEN V3.8 - 4.1: By default, reads with MAPQ < 1 and bases with BQ < 0 are ignored."
            " You can use the qc-coverage-filters-n option to specify which BQ bases and MAPQ reads to filter out.\n\n"
            "Considering only bases usable for variant calling, _i.e._ excluding:\n\n"
            "1. Clipped bases\n\n"
            "2. Bases in duplicate reads\n\n"
            "3. Reads with `MAPQ` < `min MAPQ` (default `20`)\n\n"
            "4. Bases with `BQ` < `min BQ` (default `10`)\n\n"
            "5. Reads with `MAPQ` = `0` (multimappers)\n\n"
            "6. Overlapping mates are double-counted\n\n"
            "Each _coverage_metrics.csv file may have an associated _overall_mean_cov.csv file. "
            "The latter contains the 'Average alignment coverage over &#60;source file&#62;' metric. "
            " Information about &#60;source file&#62;s can be found in the section's description"
            " or in this drop-list below if the produced text is long."
            " If input directory does not contain _overall_mean_cov files, then "
            "\"No 'coverage bed/target bed/wgs' source file found\" is printed."
        )
        for phenotype in plots:
            # Make section only if headers are not empty.
            if plots[phenotype]["headers"]:
                region_text = ""
                while regions[phenotype]:
                    region_text += improve_region(regions[phenotype].pop()) + ", "
                region_text = region_text[:-2] + ". "
                description = (
                    "Coverage metrics over "
                    + region_text
                    + bed_texts[phenotype]["description"]
                    + "\n\nPress the `Help` button for details."
                )
                anchor = "dragen-cov-metrics-own-section-" + re.sub(r"(\s|-|\.|_)+", "-", phenotype)
                plots[phenotype]["config"]["id"] = f"{anchor}-table"
                sections.append(
                    {
                        "name": make_section_name(phenotype),
                        "anchor": anchor,
                        "helptext": helptext + bed_texts[phenotype]["helptext"],
                        "description": description,
                        "plot": table.plot(
                            plots[phenotype]["data"],
                            clean_headers(order_headers(plots[phenotype]["headers"])),
                            plots[phenotype]["config"],
                        ),
                    }
                )
        return sections

    # Return the closures.
    return make_general_stats, make_own_coverage_sections


make_general_stats, make_cov_sections = create_table_handlers()


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-prefix>.<coverage-region-prefix>_coverage_metrics<arbitrary-suffix>.csv

- invalid_file_lines, which do not conform to:
    COVERAGE SUMMARY,,<metric>,<value1>  or
    COVERAGE SUMMARY,,<metric>,<value1>,<value2>

debug/info:
- unknown_metrics are those, which could not be recognized by the code.
  Their headers are incomplete, hence uninformative and ugly table's columns.
- unusual_values are those except for int/float/NA.
'''
log_data = {
    "invalid_file_names": defaultdict(list),
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_metrics": [],
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
}


def construct_coverage_parser():
    """Isolation for the parsing codeblock. Returns the closure coverage_metrics_parser."""

    '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    The DRAGEN v3.7 standard does not provide any info about consistency of coverage metrics.
    In practice two metrics can be inconsistent:
        PCT of region with coverage [ix, inf)
        PCT of region with coverage [ix, jx)

    Presumably some part of the DRAGEN tries to improve the appearance of content in csv files.
    The two mentioned metrics may have a different number of spaces inserted within them.
    Amount and insertion places seem to depend on the highest number of digits in ix and jx.
    So, for example:
        if max ix = jx = 100x   [100x: inf), [ 50x:100x)   then [  0x: inf),  [ 10x: inf)
                                                            and [  0x:  1x),  [ 10x: 15x)
    but if max ix = jx = 1500x  [1500x: inf),[1000x:1500x) then [   0x: inf), [  10x: inf)
                                                            and [   0x:   1x),[  10x:  15x)

    Using such versions of the same metric as keys would result in two columns instead of one.
    There is a possibility that an input folder contains a mixed data with different versions.

    Assumption: each metric may be inconsistent due to unknown reasons.
    Solution:   modify each metric.

    To not exaggerate the problem, only the most obvious possible issues are solved.
    If new metric-specific peculiarities arise, they can be added to the make_consistent_metric.
    '''

    PCT_RGX = re.compile("^(PCT of .+ with coverage )(.+)", re.IGNORECASE)

    def make_consistent_metric(metric):
        """Tries to fix consistency issues that may arise in coverage metrics data."""
        metric = re.sub(r"\s+", " ", metric).strip()

        pct_case = PCT_RGX.search(metric)
        if pct_case:
            metric = pct_case.group(1) + pct_case.group(2).replace(" ", "")

        metric_id = metric.lower()

        # metric will be stored in .txt output.
        # metric_id will be used for internal storage.
        return metric, metric_id

    def modify_arbitrary_suffix(suffix):
        """Support function to handle <arbitrary-suffix>."""
        # If you wish to modify this function, then you also have to adapt overall_mean_cov.py

        # Return empty string for 'tumor'.
        if re.search("tumor", suffix, re.IGNORECASE):
            return ""
        else:
            return suffix

    FILE_RGX = re.compile(r"(.+)\.(.+)_coverage_metrics(.*)\.csv$")
    LINE_RGX = re.compile(r"^COVERAGE SUMMARY,,([^,]+),([^,]+),?([^,]*)$")

    def coverage_metrics_parser(file_handler):
        """Parser for coverage metrics csv files.
        Input:  file_handler with necessary info - file name/content/root
        Output: {"success": False} if file name test failed. Otherwise:
                {"success": False/True,
                 "sample_name": <output-prefix> + <arbitrary-suffix>,
                 "phenotype": <coverage-region-prefix>,
                 "data": {
                    "data": extracted metrics with values
                    "region": extracted region (eg genome)
                 },
                 "metric_IDs_with_original_names": {
                    metric_ID: original_metric_name
                 },
                }, with success False if all file's lines are invalid.
        """
        """ File name is checked below.
        Because there is no guarantee that it would be as in the official standard:
        <output-prefix>.<coverage-region-prefix>_coverage_metrics.csv

        Accepted structure of files:
        <output-prefix>.<coverage-region-prefix>_coverage_metrics<arbitrary-suffix>.csv

        Some possible file names:

        T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
        T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv
        T_SRR7890936_50pc.target_bed_coverage_metrics.csv
        SAMPLE_0_dragen.wgs_coverage_metrics.csv
        SAMPLE_0_dragen.qc-coverage-region-1_coverage_metrics.csv
        SAMPLE_0_dragen.qc-coverage-region-2_coverage_metrics.csv
        """
        file, root = file_handler["fn"], file_handler["root"]

        file_match = FILE_RGX.search(file)
        if not file_match:
            log_data["invalid_file_names"][root].append(file)
            return {"success": False}

        sample, phenotype, suffix = file_match.groups()
        if suffix:
            sample += modify_arbitrary_suffix(suffix)

        # Stores metric IDs with original metric strings.
        # Hopefully there are no DRAGEN version-dependent differences (eg lower case).
        metric_IDs_with_original_names = {}

        success = False
        data = {}
        for line in file_handler["f"].splitlines():
            # Check the general structure. Maybe a new section will be added in the future.
            line_match = LINE_RGX.search(line)

            # If line is fine then extract the necessary fields.
            if line_match:
                success = True
                metric, value1, value2 = line_match.groups()

            # Otherwise check if line is empty. If not then report it. Go to the next line.
            else:
                if not re.search(r"^\s*$", line):
                    log_data["invalid_file_lines"][root][file].append(line)
                continue

            consistent_metric, consistent_metric_id = make_consistent_metric(metric)
            metric_IDs_with_original_names[consistent_metric_id] = consistent_metric

            """ Check the extracted values.
            These are int/float and in some cases NA for now. All other values will be reported.
            The type conversion is only performed to int and float. Other types of data are left
            unchanged as strings. It has consequences for the "modify" configuration, which are
            described in the comment for the METRICS.
            No additional try blocks are required to check other simple types of values(eg None).
            Just adapt the regexes.
            """
            try:
                value1 = int(value1)
            except ValueError:
                try:
                    value1 = float(value1)
                except ValueError:
                    if not re.search(r"^NA$", value1.strip(), re.IGNORECASE):
                        log_data["unusual_values"][root][file][metric] = value1

            data[consistent_metric_id] = value1

            if value2:
                try:
                    value2 = float(value2)
                except ValueError:
                    try:
                        value2 = int(value2)
                    except ValueError:
                        if not re.search(r"^NA$", value2.strip(), re.IGNORECASE):
                            log_data["unusual_values"][root][file][metric + " (second value)"] = value2

                data[consistent_metric_id + V2] = value2

        return {
            "success": success,
            "sample_name": sample,
            "phenotype": phenotype,
            "data": {
                "data": data,
                "region": extract_coverage_region(data),
            },
            "metric_IDs_with_original_names": metric_IDs_with_original_names,
        }

    return coverage_metrics_parser  # Return the closure.


coverage_parser = construct_coverage_parser()


def create_coverage_headers_handler():
    """Isolation for all the headers-building machinery."""

    # All regexes are constructed to be as general and simple as possible to speed up the matching.
    # "make_configs" will point later to a function, which automatically sets some header's configs.
    R_B_PAT = {
        "RGX": re.compile(r"^Aligned (?P<entity>reads|bases)$", re.IGNORECASE),
        "make_configs": None,
    }
    ALN_PAT = {
        "RGX": re.compile(r"^Aligned (?P<entity>.+?) in (?P<region>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    AVG_PAT = {
        "RGX": re.compile(r"^Average (?P<entity>.+?) coverage over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    PCT_PAT = {
        "RGX": re.compile(r"^PCT of (?P<region>.+?) with coverage (?P<entity>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    UNI_PAT = {
        "RGX": re.compile(r"^Uniformity of coverage (?P<entity>.+?) over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    MED_PAT = {
        "RGX": re.compile(r"^Median (?P<entity>.+?) coverage over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    RAT_PAT = {
        "RGX": re.compile(r"(?P<entity>.+?) ratio over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    ANY_PAT = {"RGX": re.compile(r".+"), "make_configs": None}

    # The order is based on the structure of regexes and amount of patterns found in examined data.
    METRIC_PATTERNS_FOR_REGIONS = [PCT_PAT, AVG_PAT, ALN_PAT, UNI_PAT, MED_PAT, RAT_PAT]

    def extract_coverage_region(metric_IDs):
        """Extracts region (eg genome, QC coverage region, target region)
        from provided metric_IDs. In case region can not be extracted the
        "Unknown region" is returned."""

        # First try to extract from well-known metrics.
        # This block shall catch regions in most cases.
        for metric_id in metric_IDs:
            for MPAT in METRIC_PATTERNS_FOR_REGIONS:
                metric_match = MPAT["RGX"].search(metric_id)
                if metric_match:
                    return metric_match["region"]

        # Well, first step did not work. Check if it is a part of the REGIONS.
        for metric_id in metric_IDs:
            # Tries to find the well-known regions (REGIONS) in the metric.
            region = extract_region(metric_id)
            if region:
                return region

        # Previous steps did not work.
        return "Unknown region"

    # The order is based on the structure of regexes and amount of patterns found in examined data.
    METRIC_PATTERNS_FOR_HEADERS = [PCT_PAT, AVG_PAT, ALN_PAT, R_B_PAT, UNI_PAT, MED_PAT, RAT_PAT, ANY_PAT]

    def make_coverage_headers(metric_IDs):
        headers = {}
        for metric_id in metric_IDs:
            region = None
            for MPAT in METRIC_PATTERNS_FOR_HEADERS:
                metric_match = MPAT["RGX"].search(metric_id)
                if metric_match:
                    output_auto = MPAT["make_configs"](metric_match)
                    orig_metric = metric_IDs[metric_id]

                    if "region" in output_auto and output_auto["region"]:
                        region = output_auto["region"]

                    # Metric could not be recognized.
                    if "warning" in output_auto and output_auto["warning"]:
                        log_data["unknown_metrics"].append(orig_metric)

                    # First set common configs from the SINGLE_HEADER.
                    configs = {config: value for config, value in SINGLE_HEADER.items() if value is not None}
                    configs2 = {config: value for config, value in SINGLE_HEADER.items() if value is not None}
                    configs["title"], configs["description"] = orig_metric, orig_metric
                    configs2["title"], configs2["description"] = orig_metric + V2, orig_metric + V2

                    if AUTO_CONFIGS_ENABLED:
                        configs.update(output_auto["configs"])
                        if "configs2" in output_auto:
                            configs2.update(output_auto["configs2"])

                    if USER_CONFIGS_ENABLED:
                        # Special case. Concatenate with regions to allow region-specific settings.
                        if metric_id == "aligned reads" or metric_id == "aligned bases":
                            for _region in REGIONS:
                                configs_copy = configs.copy()
                                configs_copy.update(make_user_configs(metric_id, _region))
                                headers[metric_id + _region] = configs_copy

                        # Try to set user-configs for each metric, including "aligned reads/bases".
                        # The latter will be used for regions, which are not present in REGIONS.
                        configs.update(make_user_configs(metric_id, region))
                        configs2.update(make_user_configs(metric_id + V2, region))

                    headers[metric_id] = configs
                    # Always store configs for the second value.
                    # Those, which are not needed will be ignored.
                    # Simpler code in exchange for some extra space.
                    headers[metric_id + V2] = configs2
                    break

        return headers

    def get_std_configs(configs):
        """Copies the standard real/virtual configurations from configs."""
        return {config: value for config, value in configs.items() if config in SINGLE_HEADER or config in EXTRA_HEADER}

    def make_user_configs(metric, region):
        """Creates the user-defined configurations."""
        configs = {}

        if region:
            metric = metric.replace(region, "region")
            if region in REGIONS and REGIONS[region]:
                configs.update(get_std_configs(REGIONS[region]))

        pct_case = PCT_PAT["RGX"].search(metric)
        uni_case = UNI_PAT["RGX"].search(metric)
        if pct_case:
            metric = re.sub(r"\[\d+x:", "[ix:", metric)
            metric = re.sub(r":\d+x\)", ":jx)", metric)
        if uni_case:
            metric = re.sub(r"\d+", "d", metric)

        # Now set region-specific parameters for the given metric.
        if metric in METRICS:
            _configs = METRICS[metric]
            configs.update(get_std_configs(_configs))

            if region and region in _configs and _configs[region]:
                configs.update(get_std_configs(_configs[region]))

            # Try to set (ix:inf)/(ix:jx)-specific settings.
            if pct_case:
                ix_jx_match = re.search(r"\[(\d+x):(inf|\d+x)\)", pct_case["entity"])
                if ix_jx_match:
                    ix_jx = (ix_jx_match.group(1), ix_jx_match.group(2))
                    if ix_jx in _configs["extra"]:
                        _configs = _configs["extra"][ix_jx]
                        configs.update(get_std_configs(_configs))
                        if region in _configs:
                            configs.update(get_std_configs(_configs[region]))

            # Try to match the float and set specific configs.
            if uni_case:
                float_match = re.search(r"\(pct > (\d+\.\d+)\*mean\)", uni_case["entity"])
                if float_match:
                    _float = float_match.group(1)
                    if _float in _configs["extra"]:
                        _configs = _configs["extra"][_float]
                        configs.update(get_std_configs(_configs))
                        if region in _configs:
                            configs.update(get_std_configs(_configs[region]))

        return configs

    '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    The following functions are mainly used to create configurations for headers.

    Input:
    - metric_pattern_match, which stores the result of matching a metric.

    Output:
    A dictionary with keys:
    - obligatory  configs   dict with real/logical configs from the SINGLE_HEADER.
    - obligatory  configs2  if the second value is present or
                            in case when metric can not be recognized.
    - optional    region    if present in a metric and can be extracted easily.
                            Returning it is encouraged.
    - optional    warning   can be returned if metric could not be recognized.
                            It is used for logging. Returning it is encouraged.

    Please note: modify lambda/func must check for string input.
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

    ##----------------------------------------------##
    ## This segment was taken from ./dragen/utils.py
    ## as an adoptation for the previous code version.
    read_format = "{:,.1f}"
    read_count_multiplier = config.read_count_multiplier
    if read_count_multiplier == 1:
        read_format = "{:,.0f}"

    base_format = "{:,.1f}"
    base_count_multiplier = config.base_count_multiplier
    if base_count_multiplier == 1:
        base_format = "{:,.0f}"
    elif base_count_multiplier == 0.000000001:
        base_format = "{:,.2f}"
    ##----------------------------------------------##

    '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    # This comment includes some code, which may be useful for defining
    # complex "cond_formatting_rules" and "cond_formatting_colours"

    def make_hex(red, green, blue):
        """RGB to hex str."""
        red = min(abs(red), 255)
        red = ('0' if red < 16 else "") + hex(red)[2:]
        green = min(abs(green), 255)
        green = ('0' if green < 16 else "") + hex(green)[2:]
        blue = min(abs(blue), 255)
        blue = ('0' if blue < 16 else "") + hex(blue)[2:]
        return "#" + red + green + blue

    # Different distributions' functions can be used: exp, 1/exp, Gauss, gamma ...
    # It is not even necessary to implement them, just create a chunk of data
    # in advance and store it as a variable here.

    # Below is an example of a linear change from red to green on [0.0, 100.0].
    Red_to_Green_on100_colors = [
        {str(color_n): make_hex(255-color_n, color_n, 0)}
        for color_n in range(256)
    ]
    # Over 100 is red/bad.
    Red_to_Green_on100_colors.append({"bad_over_100": "#FF0000"})

    Red_to_Green_on100_rules = {
        str(color_n): [{"gt": ((color_n)/255)*100}]
        for color_n in range(256)
    }
    # Over 100 is red/bad.
    Red_to_Green_on100_rules["bad_over_100"] = [{"gt": 100}]
    """""""""""""""""""""""""""""""""""""""""""""""""""""'''

    def extract_region(metric):
        """Tries to find the well-known regions in the metric."""
        for REGION in REGIONS:
            if REGION in metric:
                return REGION
        return None

    def improve_region(region):
        """Little helper to modify the extracted lower-case region."""
        if region == "qc coverage region":
            return "QC coverage region"
        return region

    def get_Reads_Bases_configs(metric_pattern_match):
        """The following metrics are supported:
        Aligned bases, Aligned reads
        """
        if metric_pattern_match["entity"] == "bases":
            return {
                "configs": {
                    "min": 0,
                    "format": base_format,
                    "description": f"Total number ({config.base_count_desc}) of aligned bases.",
                    "title": "Aln bases",
                    "shared_key": "base_count",
                }
            }
        else:
            return {
                "configs": {
                    "min": 0,
                    "format": read_format,
                    "description": f"Total number ({config.read_count_desc}) of aligned reads.",
                    "title": "Aln reads",
                    "shared_key": "read_count",
                }
            }

    def get_Aligned_configs(metric_pattern_match):
        """The following metrics are supported:
        Aligned bases in region, Aligned reads in region
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)
        if entity == "bases":
            description = f"Number ({config.base_count_desc}) of uniquely mapped bases to " + region
            configs = {
                "min": 0,
                "format": base_format,
                "title": "Bases on target",
                "description": description + ".",
                "shared_key": "base_count",
            }
            configs2 = {
                "min": 0,
                "max": 100,
                "suffix": "%",
                "format": base_format,
                "title": "Bases on target",
                "description": description + " relative to the number of uniquely mapped bases to the genome.",
            }
        elif entity == "reads":
            description = (
                f"Number ({config.read_count_desc}) of uniquely mapped reads to "
                + region
                + ". DRAGEN V3.4 - V3.8:"
                + " When region is the target BED, this metric is equivalent to and replaces"
                + " Capture Specificity based on target region."
                + " DRAGEN V3.9 - V4.1: Only reads with MAPQ >= 1 are included."
                + " Secondary and supplementary alignments are ignored."
            )
            configs = {
                "min": 0,
                "format": read_format,
                "title": "Reads on target",
                "description": description,
                "shared_key": "read_count",
            }
            configs2 = {
                "min": 0,
                "max": 100,
                "suffix": "%",
                "format": read_format,
                "title": "Reads on target",
                "description": "Number of uniquely mapped reads to "
                + region
                + " relative to the number of uniquely mapped reads to the genome.",
            }
        # Else unrecognized.
        else:
            metric = metric_pattern_match.string
            return {
                "configs": {
                    "title": "Aln " + entity,
                    "description": metric,
                    "scale": "Purples",
                },
                "configs2": {
                    "title": "Aln " + entity + V2,
                    "description": metric + ". Second value.",
                    "scale": "Reds",
                },
                "warning": True,
                "region": extract_region(metric),
            }

        return {"configs": configs, "configs2": configs2, "region": REGION}

    def get_Average_configs(metric_pattern_match):
        """The following metrics are supported:
        Average alignment coverage over region
        Average chromosome X/Y coverage over region
        Average mitochondrial coverage over region
        Average autosomal coverage over region
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)

        configs = {
            "suffix": " x",
            "min": 0,
        }
        if entity == "alignment":
            configs["title"] = "Depth"
            configs["description"] = (
                "Coverage depth over "
                + region
                + ": number of uniquely mapped bases to "
                + region
                + " divided by the number of sites in "
                + region
                + "."
            )
        elif entity == "autosomal":
            configs["title"] = "AvgAutCov"
            configs["description"] = (
                "Total number of bases that aligned to the autosomal loci in "
                + region
                + " divided by the total number of loci in the autosomal loci in "
                + region
                + ". If there is no autosome in the reference genome, or the "
                + region
                + " does not intersect autosomes, this metric shows as NA."
            )
        elif entity == "mitochondrial":
            configs["title"] = "AvgMitoCov"
            configs["description"] = (
                "Total number of bases that aligned to the intersection of the mitochondrial chromosome with "
                + region
                + " divided by the total number of loci in the intersection of the mitochondrial chromosome with "
                + region
                + ". If there is no mitochondrial chromosome in the reference genome or the "
                + region
                + " does not intersect mitochondrial chromosome, this metric shows as NA."
            )
        elif entity == "chr x":
            configs["title"] = "AvgXcov"
            configs["description"] = (
                "Total number of bases that aligned to the intersection of chromosome X with "
                + region
                + " divided by the total number of loci in the intersection of chromosome X with "
                + region
                + ". If there is no chromosome X in the reference genome or the "
                + region
                + " does not intersect chromosome X this metric shows as NA."
            )
        elif entity == "chr y":
            configs["title"] = "AvgYcov"
            configs["description"] = (
                "Total number of bases that aligned to the intersection of chromosome Y with "
                + region
                + " divided by the total number of loci in the intersection of chromosome Y with "
                + region
                + ". If there is no chromosome Y in the reference genome or the "
                + region
                + " does not intersect chromosome Y this metric shows as NA."
            )
        # Else unrecognized.
        else:
            return {
                "configs": {
                    "title": "Average " + entity,
                    "description": metric_pattern_match.string,
                    "scale": "Purples",
                },
                "configs2": {
                    "title": "Average " + entity + V2,
                    "description": metric_pattern_match.string + ". Second value.",
                    "scale": "Reds",
                },
                "warning": True,
                "region": REGION,
            }
        return {"configs": configs, "region": REGION}

    # Special regex to match (PCT > d.d*mean) substring.
    # Placed here because it may be used more than once.
    RGX_FLOAT = re.compile(r"\(pct\s*>\s*(\d+\.\d+)\*mean\)")

    def get_Uniformity_configs(metric_pattern_match):
        """The following metrics are supported:
        Uniformity of coverage (PCT > d.d*mean) over region
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)

        entity_match = RGX_FLOAT.search(entity)
        if entity_match:
            multiplier = entity_match.group(1)
            percent = str(float(multiplier) * 100) + "%"
            configs = {
                "suffix": "%",
                "min": 0,
                "max": 100,
                "title": "Uniformity(>" + multiplier + "&#215;mean)",
                "description": "Percentage of sites with coverage greater than "
                + percent
                + " of the mean coverage in "
                + region
                + ".",
            }
            if multiplier == "0.2":
                configs["description"] += " Demonstrates the uniformity of coverage, the higher the better."

            return {"configs": configs, "region": REGION}

        # Else unrecognized.
        return {
            "configs": {
                "title": "Uniformity " + entity,
                "description": metric_pattern_match.string,
                "scale": "Purples",
            },
            "configs2": {
                "title": "Uniformity " + entity + V2,
                "description": metric_pattern_match.string + ". Second value.",
                "scale": "Reds",
            },
            "warning": True,
            "region": REGION,
        }

    # Special regex to match the [ix, inf) and [ix, jx)
    # Placed here because it will be used more than once.
    IX_JX = re.compile(r"\[(\d+)x:(inf|\d+x)\)")

    def get_PCT_configs(metric_pattern_match):
        """The following metrics are supported:
        PCT of region with coverage [ix, inf)
        PCT of region with coverage [ix, jx)
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)

        entity_match = IX_JX.search(entity)
        if entity_match:
            IX = entity_match.group(1)
            JX = entity_match.group(2)
            description = "Percentage of sites in " + region
            if JX == "inf":
                title = "≥" + IX + "x"

                if IX == "0":
                    description += " with any coverage."
                else:
                    description += " with at least " + IX + "x coverage."
            else:
                JX = JX[:-1]  # Get rid of 'x' substring.
                title = entity
                if IX == "0" and JX == "1":
                    description += " with no coverage."
                else:
                    description += " with at least " + IX + "x but less than " + JX + "x coverage."

            return {
                "configs": {
                    "min": 0,
                    "max": 100,
                    "suffix": "%",
                    "title": title,
                    "description": description,
                },
                "region": REGION,
            }
        # Else unrecognized, eg. [ix, jx]
        return {
            "configs": {
                "title": "PCT " + entity,
                "description": metric_pattern_match.string,
                "scale": "Purples",
            },
            "configs2": {
                "title": "PCT " + entity + V2,
                "description": metric_pattern_match.string + ". Second value.",
                "scale": "Reds",
            },
            "warning": True,
            "region": REGION,
        }

    def get_Median_configs(metric_pattern_match):
        """The following metrics are supported:
        Median autosomal coverage over region
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)

        if entity == "autosomal":
            return {
                "configs": {
                    "suffix": " x",
                    "min": 0,
                    "title": "Med aut cov",
                    "description": "Median alignment coverage over the autosomal loci in "
                    + region
                    + ". If there is no autosome in the reference genome or the "
                    + region
                    + " does not intersect autosomes, this metric shows as NA.",
                },
                "region": REGION,
            }
        # Else unrecognized.
        return {
            "configs": {
                "title": "Med " + entity,
                "description": metric_pattern_match.string,
                "scale": "Purples",
            },
            "configs2": {
                "title": "Med " + entity + V2,
                "description": metric_pattern_match.string + ". Second value.",
                "scale": "Reds",
            },
            "warning": True,
            "region": REGION,
        }

    def get_ratio_over_configs(metric_pattern_match):
        """The following metrics are supported:
        Mean/Median autosomal coverage ratio over region
        XAvgCov/YAvgCov ratio over region
        XAvgCov/AutosomalAvgCov ratio over region
        YAvgCov/AutosomalAvgCov ratio over region
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)

        if entity == "mean/median autosomal coverage":
            return {
                "configs": {
                    "suffix": " x",
                    "title": "Mean/Med autosomal coverage",
                    "min": 0,
                    "description": "Mean autosomal coverage in "
                    + region
                    + " divided by the median autosomal coverage in "
                    + region
                    + ". If there is no autosome in the reference genome or the "
                    + region
                    + " does not intersect autosomes, this metric shows as NA.",
                },
                "region": REGION,
            }

        elif entity == "xavgcov/yavgcov":
            return {
                "configs": {
                    "min": 0,
                    "title": "XAvgCov/YAvgCov",
                    "description": "Average chromosome X alignment coverage in "
                    + region
                    + " divided by the average chromosome Y alignment coverage in "
                    + region
                    + ". If there is no chromosome X or chromosome Y in the reference genome or the "
                    + region
                    + " does not intersect chromosome X or Y, this metric shows as NA.",
                },
                "region": REGION,
            }

        elif entity == "xavgcov/autosomalavgcov":
            return {
                "configs": {
                    "min": 0,
                    "title": "XAvgCov/AutAvgCov",
                    "description": "Average chromosome X alignment coverage in "
                    + region
                    + " divided by the average autosomal coverage in "
                    + region
                    + ".",
                },
                "region": REGION,
            }

        elif entity == "yavgcov/autosomalavgcov":
            return {
                "configs": {
                    "min": 0,
                    "title": "XAvgCov/AutAvgCov",
                    "description": "Average chromosome Y alignment coverage in "
                    + region
                    + " divided by the average autosomal coverage in "
                    + region
                    + ".",
                },
                "region": REGION,
            }

        # Else unrecognized.
        else:
            return {
                "configs": {
                    "title": entity,
                    "description": metric_pattern_match.string,
                    "scale": "Purples",
                },
                "configs2": {
                    "title": entity + V2,
                    "description": metric_pattern_match.string + ". Second value.",
                    "scale": "Reds",
                },
                "warning": True,
                "region": REGION,
            }

    def get_any_configs(metric_pattern_match):
        """Set a bare minimum for other Metrics."""
        metric = metric_pattern_match.string
        return {
            "configs": {
                "title": metric,
                "description": metric,
                "scale": "Purples",
            },
            "configs2": {
                "title": metric + V2,
                "description": metric + ". Second value.",
                "scale": "Reds",
            },
            "warning": True,
            "region": extract_region(metric),
        }

    R_B_PAT["make_configs"] = get_Reads_Bases_configs
    ALN_PAT["make_configs"] = get_Aligned_configs
    AVG_PAT["make_configs"] = get_Average_configs
    PCT_PAT["make_configs"] = get_PCT_configs
    UNI_PAT["make_configs"] = get_Uniformity_configs
    MED_PAT["make_configs"] = get_Median_configs
    RAT_PAT["make_configs"] = get_ratio_over_configs
    ANY_PAT["make_configs"] = get_any_configs

    # Return the closures.
    return make_coverage_headers, extract_coverage_region


make_coverage_headers, extract_coverage_region = create_coverage_headers_handler()
