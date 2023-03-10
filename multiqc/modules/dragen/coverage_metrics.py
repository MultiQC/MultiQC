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

import re
import logging
from collections import defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file

from .utils import order_headers, clean_headers

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN coverage"


'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The following code section provides a clean and simple tool for setting configurations:
https://github.com/ewels/MultiQC/blob/master/docs/plots.md#creating-a-table

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
    "minRange": None,  # Minimum range for automatic bar
    "scale": "GnBu",  # Colour scale for colour coding. False to disable.
    "bgcols": None,  # Dict with values: background colours for categorical data.
    "colour": "15, 150, 255",  # Colour for column grouping
    "suffix": "",  # Suffix for value (eg. "%")
    "format": "{:,.2f}",  # Value format string - MultiQC default 1 decimal places
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values, special case, see below
    "hidden": True,  # Set to True to hide the column in the general table on page load.
    "hidden_own": False,  # For non-general plots in own coverage sections.
    # "exclude": True,         # Exclude all headers from the general html table.
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
    "hidden_own": False,  # For non-general plots in the coverage section.
    "exclude": False,  # True to exclude metric from the general html table.
    "exclude_own": False,  # True to exclude from own coverage section html table.
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
# The desired regions can be added:
# TRGT = "target region"
BED = "target bed"  # May be incorrect. No real data was available.

"""
The REGIONS holds the variables. You can add common region's settings here.
These configs have higher priority than those in the SINGLE_HEADER.
They also overwrite the automatically created ones.
REGIONS is also used to detect actual regions from real data, so
it is a good idea to include all known regions in the dictionary.
"""
REGIONS = {
    WGS: {
        # "hidden": False,
        "scale": "Blues",
    },
    QC: {
        # "exclude": True,
        "scale": "Greens",
    },
    BED: {
        "scale": "Reds",
    },
}

'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The METRICS is the container for "high level" representations of the standard metrics.
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
  * Single backslashes must be escaped.
  * If it still does not work, then you might need to take a look
    at the "make_metric_id" and maybe at the "make_user_configs".

- Some "title"s below contain the "&nbsp;" html entity to widen the column's title.
  Otherwise long values can overlap with adjacent values in right columns/rows below.
  Basically it is just a dirty little trick, which can be useful in rare cases.
  These entities will be deleted for the general table, since its columns' titles
  contain extra region suffixes, which make the columns wide enough.

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
  to none. "Configure Columns" is not included in table's div block.

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

- The problem with overlapping adjacent columns/rows, when column's
  title block is shorter than single value blocks below.
  There is a solution, which can automatically improve the appearance
  of html tables. After data/headers have been collected, a special
  function converts the column's value to string, concatenates with
  "suffix" if presented and appends at least abs(maxlen(values) - len(title))
  "&nbsp;" entities. Simple, but does not solve the issue in general.
  Additional JS script could analyze tables' data and attributes(eg font)
  to decide how much extra padding space is needed to increase table's
  readability further to some extent.


Final notes:

- Most configs below are for demonstration purposes only.
  You can adapt them according to your personal taste/needs.

- You have the final word in defining settings.
  The code trusts you and does not check the chosen values.
"""

# These might be useful for checking results.
AUTO_CONFIGS_ENABLED = True
USER_CONFIGS_ENABLED = True

V2 = "(second value)"  # Used to create a unique ID if a metric has two values.

""" GTQ is the "greater than or equal to" sign.
Other types can be found in the "Extra Info" comment above.
GTQ is used by setting automatic titles for the metric:
PCT of region with coverage [ix: inf)
"""
GTQ = "&ge;"

METRICS = {
    "average chr x coverage over region": {
        "order_priority": 1.0,
        "title": "X cov",
        "suffix": " x",
        "colour": "255, 255, 0",
        WGS: {
            "hidden": False,
        },
    },
    "average chr y coverage over region": {
        "order_priority": 1.1,
        "title": "Y cov",
        "suffix": " x",
        "scale": "YlOrRd",
        "colour": "255, 255, 0",
        WGS: {
            "hidden": False,
        },
    },
    "xavgcov/yavgcov ratio over region": {
        "order_priority": 1.2,
        "title": "XAvgCov/YAvgCov",
        "hidden_own": True,
    },
    "xavgcov/autosomalavgcov ratio over region": {
        "order_priority": 1.3,
        "title": "XAvgCov/AutAvgCov",
        "hidden_own": True,
    },
    "yavgcov/autosomalavgcov ratio over region": {
        "order_priority": 1.4,
        "title": "YAvgCov/AutAvgCov",
        "hidden_own": True,
    },
    "average alignment coverage over region": {
        "order_priority": 3,
        "title": "Depth",
        "scale": "BrBG",
        "colour": "0, 255, 255",
    },
    "average mitochondrial coverage over region": {
        "order_priority": 4,
        "title": "MT cov&nbsp;&nbsp;",
        "scale": "Reds",
        "suffix": " x",
        "colour": "255, 0, 255",
    },
    "average autosomal coverage over region": {
        "order_priority": 6,
        "title": "Mean aut cov",
        "suffix": " x",
        "scale": "Reds",
        "bgcols": {"NA": "#FF00FF"},
        WGS: {
            "description": "Average autosomal coverage over genome. Calculated "
            + "as the number of bases that aligned to the autosomal loci in genome"
            + " divided by the total number of loci in the autosomal loci in "
            + "genome. If there are no autosomes in the reference genome, or "
            + "the region does not intersect autosomes, this metric shows as NA.",
        },
    },
    "median autosomal coverage over region": {
        "order_priority": 7,
        "title": "Med aut cov",
        "suffix": " x",
    },
    "mean/median autosomal coverage ratio over region": {
        "order_priority": 7.1,
        "title": "Mean/med aut cov",
        "suffix": " x",
        "hidden_own": True,
    },
    # "aligned bases" has no support for region-specificity.
    # So, you can not set WGS/QC/BED-specific settings here.
    # Well technically you can but it uses the last found region.
    "aligned bases": {
        "order_priority": 8,
        "title": "Aln bases",
        "scale": "RdYlGn",
        "colour": "0, 0, 255",
    },
    "aligned bases in region": {
        "order_priority": 9,
        "title": "Bases on target",
        "scale": "Reds",
        "colour": "0, 0, 255",
        WGS: {
            "scale": "Purples",
        },
        QC: {
            "scale": "Oranges",
        },
    },
    "aligned bases in region"
    + V2: {
        "order_priority": 10,
        "title": "Bases on trg pct",
        "max": 100,
        "suffix": " %",
        "bgcols": {
            "NA": "#00FFFF",
        },
    },
    # The same as "aligned bases".
    "aligned reads": {
        "order_priority": 11,
        "title": "Aln reads",
    },
    "aligned reads in region": {
        "order_priority": 12,
        "title": "Reads on target",
        "scale": "RdGy",
        "colour": "255, 0, 0",
    },
    "aligned reads in region"
    + V2: {
        "order_priority": 13,
        "title": "Reads on trg pct",
        "max": 100,
        "suffix": " %",
        "scale": "RdGy",
        "colour": "255, 0, 0",
    },
    # This metric can handle any float. "0.0" is just a technical trick.
    # Only the ">" operator and the "*mean" are supported.
    # Title and description shall not be set here, they are automated.
    # But if you want to handle a certain number manually, then you can
    # add a number-string key with desired configs in the "extra".
    "uniformity of coverage (pct > 0.0*mean) over region": {
        "suffix": " %",
        "colour": "55, 255, 55",
        "extra": {
            # "Uniformity of coverage (PCT > 0.2*mean) over region" metric.
            "0.2": {
                "scale": "PiYG",
                "order_priority": 5.2,
                # "Uniformity of coverage (PCT > 0.2*mean) over QC coverage region" metric.
                QC: {
                    "description": "Percentage of sites with coverage greater than "
                    + "20% of the mean coverage in QC coverage region. Demonstrates"
                    + " the uniformity of coverage, the higher the better.",
                },
                # "Uniformity of coverage (PCT > 0.2*mean) over genome" metric.
                WGS: {
                    "description": "Percentage of sites with coverage greater "
                    + "than 20% of the mean coverage in genome. Demonstrates"
                    + " the uniformity of coverage, the higher the better.",
                },
            },
            # "Uniformity of coverage (PCT > 0.4*mean) over region" metric.
            "0.4": {
                "order_priority": 5.4,
                QC: {
                    "description": "Percentage of sites with coverage greater "
                    + "than 40% of the mean coverage in QC coverage region.",
                },
            },
        },
    },
    # The following two metrics are the most tedious. The creation of the title
    # and description are automated. If you want to target some specific case,
    # then you can specify a certain combination of (ix, inf) or (ix, jx).
    # Just write a group in the "extra", you can consider it as a separate
    # metric, so you can also add general or region-specific configs there.
    # "0x" are just a tecnical trick. Any natural number is meant instead.
    # "PCT of region with coverage [ix, inf)" metrics:
    "pct of region with coverage [0x:inf)": {
        "max": 100,
        "scale": "Purples",
        "suffix": " %",
        "colour": "255, 50, 25",
        WGS: {
            # "scale": "Reds",
        },
        QC: {
            # "scale": "Oranges",
        },
        "extra": {
            # "PCT of region with coverage [0x: inf)" metric.
            ("0x", "inf"): {
                "title": GTQ + "0x" + 8 * "&nbsp;",
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
                BED: {
                    "description": "Percentage of sites in target bed with any coverage.",
                },
            },
            # "PCT of region with coverage [1x: inf)" metric.
            ("1x", "inf"): {
                # "title": "Some title",
                # "description": "Some description", ...
            },
            # "PCT of region with coverage [100x: inf)" metric.
            ("100x", "inf"): {
                # "title": "Another title",
                # "description": "Some description", ...
            },
        },
    },
    # "PCT of region with coverage [ix, jx)" metrics:
    "pct of region with coverage [0x:0x)": {
        "max": 100,
        "scale": "Reds",
        "suffix": " %",
        "hidden_own": True,
        # "exclude": True,
        # "exclude_own": True,
        "extra": {
            ("0x", "1x"): {
                "title": "0x&nbsp;&nbsp;&nbsp;",
                WGS: {
                    "description": "Percentage of sites in genome with no coverage.",
                },
                QC: {
                    "description": "Percentage of sites in QC coverage region with no coverage.",
                },
                BED: {
                    "description": "Percentage of sites in target bed with no coverage.",
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
    "sortRows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}
# Below are region-specific configs with higher priority.
REGION_TABLE_CONFIG = {
    WGS: {
        "table_title": "WGS",
        # "no_beeswarm": True,
        # "save_file": True,
        # "raw_data_fn": "Genome_raw_file",
    },
    QC: {
        "table_title": "QC coverage region",
    },
    BED: {
        "table_title": "Target bed",
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

    def add_coverage_metrics(self):
        """The main function of the dragen coverage metrics module.
        The public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns either a view object containing
        a list of sample names, or an empty set if data could not be collected or
        became empty after ignoring samples."""

        cov_data = defaultdict(dict)
        cov_headers = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones,
        # and to extract data from _overall_mean_cov.csv files.
        all_samples = defaultdict(lambda: defaultdict(list))

        for file in self.find_log_files("dragen/coverage_metrics"):
            out = cov_parser(file, cov_headers)
            if out["success"]:
                self.add_data_source(file, section="stats")

                original_sample = out["sample_name"]
                cleaned_sample = self.clean_s_name(original_sample, file)
                phenotype = out["phenotype"]

                all_samples[cleaned_sample][phenotype].append((original_sample, file))

                # Add/overwrite the sample.
                cov_data[cleaned_sample][phenotype] = out["data"]

        '''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Coverage data is collected now. The headers are also prepared.
        If other files (eg json) are needed, then these can be found
        with find_log_files, processed and used afterwards.
        See afterqc.py or supernova.py for example.
        CSS/JS can be imported, see fastqc.py for details.
        '''
        cov_data = self.ignore_samples(cov_data)

        if not cov_data:
            return set()

        # Check samples for duplicates.
        check_duplicate_samples(all_samples)

        # Write data into the general table.
        gen_data, gen_headers = make_general_stats(cov_data, cov_headers)
        self.general_stats_addcols(gen_data, gen_headers, namespace=NAMESPACE)

        # Write data to file.
        out_data = make_data_for_txt_report(cov_data)
        self.write_data_file(out_data, "dragen_cov_metrics")

        # Extract coverage bed/target bed/wgs from _overall_mean_cov.csv files.
        # And prepare <coverage region prefix>-specific texts.
        bed_texts = make_bed_texts(self.overall_mean_cov_data, all_samples)
        coverage_sections = make_cov_sections(cov_data, cov_headers, bed_texts)

        # Special closure for reporting found info/warnings/errors,
        # which were collected while calling the cov_parser.
        # You can disable it anytime, if it is not wanted.
        make_parsing_log_report()

        for cov_section in coverage_sections:
            self.add_section(
                name=cov_section["name"],
                anchor=cov_section["anchor"],
                description=cov_section["description"],
                helptext=cov_section["helptext"],
                plot=cov_section["plot"],
            )
        return cov_data.keys()


def check_duplicate_samples(sample_names):
    """Check samples for duplicate names. Warn about found ones."""
    message = ""
    line1 = "\n  {}.{} was built from the following samples:"
    line2 = "\n    {} in {}"
    for s_name in sample_names:
        for phenotype in sample_names[s_name]:
            if len(sample_names[s_name][phenotype]) > 1:
                message += line1.format(s_name, phenotype)
                for orig_s_name, file in sample_names[s_name][phenotype]:
                    message += line2.format(file["fn"], file["root"])
    if message:
        message = (
            "\nDuplicate sample names were found. The last one overwrites previous data."
            + message
        )
        log.warning(message)


def make_data_for_txt_report(coverage_data):
    """Prepare data for the text report."""
    data = {}
    for sample in coverage_data:
        for phenotype in coverage_data[sample]:
            ID = sample + phenotype
            data[ID] = coverage_data[sample][phenotype]["metrics_and_values"]
    return data


def make_bed_texts(overall_mean, sample_names):
    """Matches _overall_mean_cov.csv to the corresponding _coverage_metrics.csv
    Extracts coverage bed/target bed/wgs/file names and creates a text for each
    section (phenotype)."""

    # Each sample.phenotype can have at most 1 corresponding overall_mean_cov.csv file.
    # If it has it, then append the sample to the sources_matched. If not, then
    # append to the sources_not_matched, in order to properly deliver info to users.
    sources_matched = defaultdict(list)
    sources_not_matched = defaultdict(list)
    phenotypes = set()  # Stores all found phenotypes.

    for sample in sample_names:
        for phenotype in sample_names[sample]:
            phenotypes.add(phenotype)
            # Extract the last tuple from list with duplicates.
            orig_sample, file_handler = sample_names[sample][phenotype][-1]
            root = file_handler["root"]
            # Check if that file is present in the overall_mean_cov data.
            # Presumeably only in rare cases no data can be collected.
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

    texts_for_sections = {}
    # Now go through the collected phenotypes.
    for phenotype in phenotypes:
        # Does at least 1 sample has a source file?
        if phenotype in sources_matched:
            # There are probably not too many sources according to examined data.
            # Restructure first, so each source references list with samples.
            bed_sources = defaultdict(list)
            for sample, source in sources_matched[phenotype]:
                bed_sources[source].append(sample)

            text = ""
            # Maybe not all samples were matched.
            if phenotype in sources_not_matched:
                # First tell about matched ones.
                for source in bed_sources:
                    text += (
                        "The following samples are based on the "
                        + extract_source(source)
                        + ": "
                    )
                    for sample in bed_sources[source]:
                        text += sample + ", "

                    # Get rid of the last ", " and append a ".\n\n"
                    text = text[:-2] + ".\n\n"

                # Also tell about unmatched samples.
                text += "The following samples do not have associated source files: "
                for sample in sources_not_matched[phenotype]:
                    text += sample + ", "

                # Get rid of the last ", " and append a "."
                text = text[:-2] + "."

            # Each sample has some corresponding overall_mean_cov file.
            else:
                # Just 1 source file for all samples?
                if len(bed_sources) == 1:
                    text = (
                        "All samples are based on the "
                        + extract_source(list(bed_sources)[0])
                        + "."
                    )
                # At least 2 source files are present.
                else:
                    for source in bed_sources:
                        text += (
                            "The following samples are based on the "
                            + extract_source(source)
                            + ": "
                        )
                        for sample in bed_sources[source]:
                            text += sample + ", "
                        # Get rid of the last ", " and append ".\n\n"
                        text = text[:-2] + ".\n\n"
                    text = text[:-2]  # Get rid of the last "\n\n"

        # Otherwise not a single sample has a source from overall_mean_cov.csv.
        else:
            # Try to make the text more specific.
            if re.search("^wgs$", phenotype, re.IGNORECASE):
                text = "No wgs source file found."
            else:
                text = "No 'coverage bed/target bed' source files found."

        texts_for_sections[phenotype] = text

    return texts_for_sections


def create_table_handlers():
    """Codeblock for handling tables.
    Output 2 closures:
    * make_general_stats
    * make_own_coverage_sections"""

    # Regexes for the phenotypes:
    # <coverage region prefix> + <arbitrary suffix>, where
    # the former is the standard DRAGEN region prefix
    # and the latter is some string suffix (empty in most cases).
    RGX_WGS_PHENOTYPE = re.compile("^wgs(.*)", re.IGNORECASE)
    RGX_QC_PHENOTYPE = re.compile("^qc-coverage-region-(.+)", re.IGNORECASE)
    RGX_BED_PHENOTYPE = re.compile("^target_bed(.*)", re.IGNORECASE)

    def improve_gen_phenotype(phenotype):
        """Improve phenotype's appearance for columns' titles in the general table."""

        pheno_match = RGX_QC_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "QC-" + pheno_match.group(1)

        pheno_match = RGX_WGS_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "WGS" + pheno_match.group(1)

        pheno_match = RGX_BED_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "BED" + pheno_match.group(1)

        return phenotype

    def make_general_stats(coverage_data, headers):
        """Prepare data and headers for the general table."""

        gen_data = defaultdict(dict)
        gen_headers = {}
        for sample in coverage_data:
            for phenotype in coverage_data[sample]:
                data = coverage_data[sample][phenotype]["metrics_and_values"]
                region = coverage_data[sample][phenotype]["region"]
                for metric in data:
                    # Data and the corresponding header are included in the report,
                    # only if "exclude" is not presented or False/False-equivalent.
                    if not (
                        "exclude" in headers[metric] and headers[metric]["exclude"]
                    ):
                        # Make exclusive metric ID.
                        # Please notice that special signs (eg "]") are
                        # excluded when HTML IDs are created. So, for example:
                        # PCT of region with coverage [10x: 50x)
                        # PCT of region with coverage [10x: 50x]
                        # will both reference the same HTML ID.
                        m_id = "gen table_" + phenotype + "_" + metric
                        """
                        Only the sample is used as key for gen_data.
                        Sample shall not be concatenated with the phenotype,
                        in order to maintain compatability with other
                        dragen modules, which do not contain region prefixes.
                        """
                        gen_data[sample][m_id] = data[metric]
                        gen_headers[m_id] = headers[metric].copy()
                        """
                        Some modifications are necessary to improve informativeness
                        of the general table, because several/many familiar tables
                        can be combined and inserted into it. HTML entities are also
                        deleted, because the title shall become wide enough after
                        concatenating the phenotype.
                        Check if "title" is present for safety, because
                        it can be set to None in the SINGLE_HEADER or METRICS.
                        """
                        if "title" in gen_headers[m_id]:
                            gen_headers[m_id]["title"] = (
                                re.sub("&nbsp;", "", gen_headers[m_id]["title"])
                                + " "
                                + improve_gen_phenotype(phenotype)
                            )

        return gen_data, clean_headers(order_headers(gen_headers))

    def improve_own_phenotype(phenotype):
        """Improve phenotype's appearance for non-general sections' titles."""

        pheno_match = RGX_QC_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "QC coverage region " + pheno_match.group(1)

        pheno_match = RGX_WGS_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "WGS" + pheno_match.group(1)

        pheno_match = RGX_BED_PHENOTYPE.search(phenotype)
        if pheno_match:
            return "Target BED" + pheno_match.group(1)

        # Try to make it better looking a little bit.
        phenotype = re.sub(r"-|_|\.", " ", phenotype).strip()
        return phenotype

    # Regexes for the well-known regions:
    RGX_WGS_REGION = re.compile("^genome$", re.IGNORECASE)
    RGX_QC_REGION = re.compile("^QC coverage region$", re.IGNORECASE)
    RGX_BED_REGION = re.compile("^target bed$", re.IGNORECASE)

    def improve_region(region):
        """Modify the extracted lower-case region."""

        if RGX_QC_REGION.search(region):
            return "QC coverage region"
        if RGX_WGS_REGION.search(region):
            return "genome"
        if RGX_BED_REGION.search(region):
            return "target BED"

        return region

    def make_own_coverage_sections(coverage_data, coverage_headers, bed_texts):
        """Create non-general phenotype-specific sections."""

        plots = defaultdict(lambda: defaultdict(dict))
        # There is no guarantee, that all files would have the same region.
        regions = defaultdict(set)
        for sample in coverage_data:
            for phenotype in coverage_data[sample]:
                real_data = coverage_data[sample][phenotype]["metrics_and_values"]
                region = coverage_data[sample][phenotype]["region"]
                data = {sample: {}}
                headers = {}
                for metric in real_data:
                    if not (
                        "exclude_own" in coverage_headers[metric]
                        and coverage_headers[metric]["exclude_own"]
                    ):
                        m_id = "own table_" + phenotype + "_" + metric
                        data[sample][m_id] = real_data[metric]
                        headers[m_id] = coverage_headers[metric].copy()

                        if "hidden_own" in headers[m_id]:
                            headers[m_id]["hidden"] = headers[m_id]["hidden_own"]

                plots[phenotype]["data"].update(data)
                plots[phenotype]["headers"].update(headers)
                regions[phenotype].add(region)

                plots[phenotype]["config"] = {
                    config: val
                    for config, val in TABLE_CONFIG.items()
                    if val is not None
                }
                if region in REGION_TABLE_CONFIG and REGION_TABLE_CONFIG[region]:
                    plots[phenotype]["config"].update(REGION_TABLE_CONFIG[region])

        sections = []
        helptext = (
            "The following criteria are used when calculating coverage:\n\n"
            + "* Duplicate reads and clipped bases are ignored.\n\n"
            + "* DRAGEN V3.4 - 3.7: Only reads with `MAPQ` > `min MAPQ`"
            + " and bases with `BQ` > `min BQ` are considered\n\n"
            + "* DRAGEN V3.8 - 4.1: By default, reads with MAPQ < 1"
            + " and bases with BQ < 0 are ignored."
            + " You can use the qc-coverage-filters-n option to specify which"
            + " BQ bases and MAPQ reads to filter out.\n\n"
            + "Considering only bases usable for variant calling, _i.e._ excluding:\n\n"
            + "1. Clipped bases\n\n"
            + "2. Bases in duplicate reads\n\n"
            + "3. Reads with `MAPQ` < `min MAPQ` (default `20`)\n\n"
            + "4. Bases with `BQ` < `min BQ` (default `10`)\n\n"
            + "5. Reads with `MAPQ` = `0` (multimappers)\n\n"
            + "6. Overlapping mates are double-counted"
        )
        for phenotype in plots:
            # Make section only if headers are not empty.
            if plots[phenotype]["headers"]:
                region_text = ""
                while regions[phenotype]:
                    region_text += improve_region(regions[phenotype].pop()) + ", "
                region_text = region_text[:-2] + "."
                description = (
                    "Coverage metrics over "
                    + region_text
                    + " Press the `Help` button for details.\n\n"
                    + bed_texts[phenotype]
                )
                sections.append(
                    {
                        "name": improve_own_phenotype(phenotype) + " Metrics",
                        # Spaces are replaced with hyphens to pass MultiQC lint test.
                        "anchor": "dragen-cov-metrics-own-sec-"
                        + re.sub("\s+", "-", phenotype),
                        "helptext": helptext,
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


def construct_coverage_parser():
    """Isolation for all parsing machinery.
    Returns 2 closures:
    * coverage_metrics_parser: parser for coverage data stored in csv files.
    * make_log_report: log reporter for found info/warnings/errors."""

    # All regexes are constructed to be as general and simple as possible to speed up the matching.
    # "make_configs" will point later to a function, which automatically sets some header's configs.
    ALN_PAT = {
        "RGX": re.compile("^Aligned (?P<entity>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    AVG_PAT = {
        "RGX": re.compile(
            "^Average (?P<entity>.+?) coverage over (?P<region>.+)", re.IGNORECASE
        ),
        "make_configs": None,
    }
    PCT_PAT = {
        "RGX": re.compile(
            "^PCT of (?P<region>.+?) with coverage (?P<entity>.+)", re.IGNORECASE
        ),
        "make_configs": None,
    }
    UNI_PAT = {
        "RGX": re.compile(
            "^Uniformity of coverage (?P<entity>.+?) over (?P<region>.+)", re.IGNORECASE
        ),
        "make_configs": None,
    }
    MED_PAT = {
        "RGX": re.compile(
            "^Median (?P<entity>.+?) coverage over (?P<region>.+)", re.IGNORECASE
        ),
        "make_configs": None,
    }
    RAT_PAT = {
        "RGX": re.compile("(?P<entity>.+?) ratio over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None,
    }
    ANY_PAT = {"RGX": re.compile(".+"), "make_configs": None}
    # The order is based on the structure of regexes and amount of patterns found in examined data.
    METRIC_PATTERNS = [PCT_PAT, AVG_PAT, ALN_PAT, UNI_PAT, MED_PAT, RAT_PAT, ANY_PAT]

    '''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    The log_data is the container for found info, warnings and errors.
    Collected data is stored in groups to make the output log file more readable.

    warnings:
    - invalid_file_names, which do not conform to the:
        <output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv

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
        "unknown_metrics": defaultdict(lambda: defaultdict(list)),
        "unusual_values": defaultdict(lambda: defaultdict(dict)),
    }

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
    If new metric-specific peculiarities arise, they can be added to the make_metric_id.
    '''

    def make_metric_id(metric):
        """Tries to fix consistency issues that may arise in coverage metrics data."""

        # Single backslashes are escaped.
        metric_id = re.sub("\s+", " ", metric).strip().lower()

        pct_case = PCT_PAT["RGX"].search(metric_id)
        if pct_case:
            metric_id = (
                "pct of "
                + pct_case["region"]
                + " with coverage "
                + pct_case["entity"].replace(" ", "")
            )

        return metric_id

    FILE_RGX = re.compile(r"^([^.]+)\.(.+?)_coverage_metrics(.*)\.csv$")
    LINE_RGX = re.compile("^COVERAGE SUMMARY,,([^,]+),([^,]+),?([^,]*)$")

    # Stores extracted region for each found metric.
    cov_headers_support = defaultdict(dict)

    def coverage_metrics_parser(file_handler, headers):
        """Parser for coverage metrics csv files.
        Input:  file_handler with necessary info - file name/content/root
        Output: {"success": 0} if file name test failed. Otherwise:
                {"success": 0/1, "sample_name": <output prefix>,
                 "phenotype": <coverage region prefix> + <arbitrary suffix>,
                 "region": extracted_region,
                 "data": extracted_metrics_with_values
                }, with success 0 if all file's lines are invalid.
        """
        """ File name is checked below.
        Because there is no guarantee that it would be as in the official standard:
        <output prefix>.<coverage region prefix>_coverage_metrics.csv

        Accepted structure of files:
        <output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv

        Some possible file names:

        T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
        T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv
        T_SRR7890936_50pc.target_bed_coverage_metrics.csv
        SAMPLE_0_dragen.wgs_coverage_metrics.csv
        SAMPLE_0_dragen.qc-coverage-region-1_coverage_metrics.csv
        SAMPLE_0_dragen.qc-coverage-region-2_coverage_metrics.csv
        """
        file_name, root_name = file_handler["fn"], file_handler["root"]

        file_match = FILE_RGX.search(file_name)
        if not file_match:
            log_data["invalid_file_names"][root_name].append(file_name)
            return {"success": 0}

        sample, phenotype = file_match.group(1), file_match.group(2)

        # phenotype is just concatenated with an arbitrary group, because it is non-standard.
        # phenotype is used instead of sample, because 2D html plots for table data can be
        # created only if sample names are the same.
        if file_match.group(3):
            phenotype += file_match.group(3)

        success = 0
        region = ""
        data = {}
        for line in file_handler["f"].splitlines():
            # Check the general structure. Maybe a new section will be added in the future.
            line_match = LINE_RGX.search(line)

            # If line is fine then extract the necessary fields.
            if line_match:
                metric, value1, value2 = (
                    line_match.group(1),
                    line_match.group(2),
                    line_match.group(3),
                )
                success = 1

            # Otherwise check if line is empty. If not then report it and go to the next line.
            else:
                if not re.search("^\s*$", line):
                    log_data["invalid_file_lines"][root_name][file_name].append(line)
                continue

            """ Check the extracted values.
            These are int/float and in some cases NA for now. All other values will be reported.
            The type conversion is only performed to int and float. Other types of data are left
            unchanged as strings. It has some consequences for the "modify" configuration, which
            are described in the comment for the METRICS.
            No additional try blocks are required to check other simple types of values(eg Inf).
            Just adapt the regexes.
            """
            try:
                value1 = int(value1)
            except ValueError:
                try:
                    value1 = float(value1)
                except ValueError:
                    if not re.search("^NA$", value1, re.IGNORECASE):
                        log_data["unusual_values"][root_name][file_name][
                            metric + V2
                        ] = value1

            if value2:
                try:
                    value2 = float(value2)
                except ValueError:
                    try:
                        value2 = int(value2)
                    except ValueError:
                        if not re.search("^NA$", value2, re.IGNORECASE):
                            log_data["unusual_values"][root_name][file_name][
                                metric + V2
                            ] = value2

            metric_id = make_metric_id(metric)
            if metric_id not in headers:
                for MPAT in METRIC_PATTERNS:
                    metric_match = MPAT["RGX"].search(metric_id)
                    if metric_match:
                        output_auto = MPAT["make_configs"](metric_match)

                        if "region" in output_auto and output_auto["region"]:
                            region = output_auto["region"]
                            cov_headers_support[metric_id]["region"] = region

                        # Metric could not be recognized.
                        if "warning" in output_auto and output_auto["warning"]:
                            cov_headers_support[metric_id]["warning"] = 1

                        # First set common configs from the SINGLE_HEADER.
                        configs = {
                            config: value
                            for config, value in SINGLE_HEADER.items()
                            if value is not None
                        }

                        if AUTO_CONFIGS_ENABLED:
                            auto_configs = output_auto["configs"]
                            configs.update(auto_configs)

                        if USER_CONFIGS_ENABLED:
                            user_configs = make_user_configs(metric_id, region)
                            configs.update(user_configs)

                        headers[metric_id] = configs

                        if value2:
                            configs = {
                                config: value
                                for config, value in SINGLE_HEADER.items()
                                if value is not None
                            }

                            if AUTO_CONFIGS_ENABLED:
                                if "configs2" in output_auto:
                                    auto_configs = output_auto["configs2"]
                                else:
                                    auto_configs = {
                                        "title": metric,
                                        "description": metric,
                                        "scale": "Reds",
                                    }
                                configs.update(auto_configs)

                            if USER_CONFIGS_ENABLED:
                                user_configs = make_user_configs(metric_id + V2, region)
                                configs.update(user_configs)

                            headers[metric_id + V2] = configs

                        break

            # In case a second value is found, but the header was not created.
            # Can happen if metric has a variable number of values and if the first catched
            # metric has only 1 value, but the same subsequent metric has 2 values.
            # Since it's very unlikely to happen, a very simple header is constructed.
            if value2 and (metric_id + V2 not in headers):
                headers[metric_id + V2] = {
                    "title": metric,
                    "description": metric,
                    "scale": "Reds",
                }

            if not region and "region" in cov_headers_support[metric_id]:
                region = cov_headers_support[metric_id]["region"]

            if "warning" in cov_headers_support[metric_id]:
                log_data["unknown_metrics"][root_name][file_name].append(metric)

            data[metric_id] = value1
            if value2:
                data[metric_id + V2] = value2

        if not region:
            region = "Unknown region"

        return {
            "success": success,
            "sample_name": sample,
            "phenotype": phenotype,
            "data": {"metrics_and_values": data, "region": region},
        }

    def make_log_report():
        """The only purpose of this function is to create a readable and informative log output."""

        if log_data["invalid_file_names"]:
            log_message = (
                "\n\nThe file names must conform to the following structure:\n"
                + "<output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv\n\n"
                + "The following files are not valid:\n"
            )
            for root in log_data["invalid_file_names"]:
                log_message += "  " + root + ":\n"
                for file in log_data["invalid_file_names"][root]:
                    log_message += "    " + file + "\n"

            log.warning(log_message + "\n")

        if log_data["invalid_file_lines"]:
            log_message = (
                "\n\nThe lines in files must be:\n"
                + "COVERAGE SUMMARY,,<metric>,<value1> or "
                + "COVERAGE SUMMARY,,<metric>,<value1>,<value2>\n\n"
                + "The following files contain invalid lines:\n"
            )
            for root in log_data["invalid_file_lines"]:
                log_message += "  " + root + ":\n"
                for file in log_data["invalid_file_lines"][root]:
                    log_message += "    " + file + ":\n"
                    for line in log_data["invalid_file_lines"][root][file]:
                        log_message += "      " + line + "\n"

            log.warning(log_message + "\n")

        if log_data["unknown_metrics"]:
            log_message = "\n\nThe following files contain unknown metrics:\n"

            for root in log_data["unknown_metrics"]:
                log_message += "  " + root + ":\n"
                for file in log_data["unknown_metrics"][root]:
                    log_message += "    " + file + ":\n"
                    for metric in log_data["unknown_metrics"][root][file]:
                        log_message += "      " + metric + "\n"

            log.warning(log_message + "\n")

        if log_data["unusual_values"]:
            log_message = (
                "\n\nAll metrics' values except int, float and NA are non-standard.\n"
                + "The following files contain non-standard values:\n"
            )
            for root in log_data["unusual_values"]:
                log_message += "  " + root + ":\n"
                for file in log_data["unusual_values"][root]:
                    log_message += "    " + file + ":\n"
                    for metric in log_data["unusual_values"][root][file]:
                        log_message += (
                            "      "
                            + metric
                            + " = "
                            + log_data["unusual_values"][root][file][metric]
                            + "\n"
                        )

            log.warning(log_message + "\n")

    def get_std_configs(configs_dict):
        """Copies the standard real/virtual configurations from configs_dict."""
        return {
            config: value
            for config, value in configs_dict.items()
            if config in SINGLE_HEADER or config in EXTRA_HEADER
        }

    def make_user_configs(metric, region):
        """Creates the user-defined configurations."""
        configs = {}

        # If region is not empty, in REGIONS and defined.
        if region and region in REGIONS and REGIONS[region]:
            metric = metric.replace(region, "region")
            configs.update(get_std_configs(REGIONS[region]))

        pct_case = PCT_PAT["RGX"].search(metric)
        uni_case = UNI_PAT["RGX"].search(metric)
        if pct_case or uni_case:
            metric = re.sub("\d+", "0", metric)

        # Now set region-specific parameters for the given metric.
        if metric in METRICS:
            user_configs = METRICS[metric]
            configs.update(get_std_configs(user_configs))

            if region and region in user_configs and user_configs[region]:
                configs.update(get_std_configs(user_configs[region]))

            # Try to set (ix:inf)/(ix:jx)-specific settings.
            if pct_case and "extra" in user_configs and user_configs["extra"]:
                ix_jx = re.search(r"\[(\d+x):(inf|\d+x)\)", pct_case["entity"])
                if ix_jx:
                    _group = (ix_jx.group(1), ix_jx.group(2))

                    if _group in user_configs["extra"]:
                        user_configs = user_configs["extra"][_group]
                        configs.update(get_std_configs(user_configs))
                        if region in user_configs:
                            configs.update(get_std_configs(user_configs[region]))

            # Try to match the float and set specific configs.
            if uni_case and "extra" in user_configs and user_configs["extra"]:
                entity_match = re.search(
                    r"\(PCT > (\d+\.\d+)\*mean\)", uni_case["entity"], re.IGNORECASE
                )
                if entity_match:
                    _float = entity_match.group(1)
                    if _float in user_configs["extra"]:
                        user_configs = user_configs["extra"][_float]
                        configs.update(get_std_configs(user_configs))
                        if region in user_configs:
                            configs.update(get_std_configs(user_configs[region]))

        return configs

    '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    The following functions are mainly used to create configurations for headers.

    Input:
    - metric_pattern_match, which stores the result of matching a metric.

    Output:
    A dictionary with keys:
    - obligatory  configs   dict with real/logical configs from the SINGLE_HEADER.
    - obligatory  configs2  if the second value is presented or
                            in case when metric can not be recognized.
    - optional    region    if presented in a metric and can be extracted easily.
                            Returning it is encouraged.
    - optional    warning   can be returned if metric could not be recognized.
                            The associated value is irrelevant(shall be equal to True).

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

    # Different distributions' functions can be used: exp, Poisson, Gauss, gamma ...
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
        if re.search("QC coverage region", region, re.IGNORECASE):
            return "QC coverage region"
        return region

    def get_Aligned_configs(metric_pattern_match):
        """The following metrics are supported:
        Aligned bases, Aligned reads
        Aligned bases in region, Aligned reads in region
        """
        entity = metric_pattern_match["entity"]

        if entity == "bases" or entity == "reads":
            if entity == "bases":
                return {
                    "configs": {
                        "min": 0,
                        "format": base_format,
                        "description": "Total number ({}) of aligned bases.".format(
                            config.base_count_desc
                        ),
                        "title": config.base_count_prefix + " Aln bases",
                        "modify": lambda x: x
                        if isinstance(x, str)
                        else x * base_count_multiplier,
                    }
                }
            else:
                return {
                    "configs": {
                        "min": 0,
                        "format": read_format,
                        "description": "Total number ({}) of aligned reads.".format(
                            config.read_count_desc
                        ),
                        "title": config.read_count_prefix + " Aln reads",
                        "modify": lambda x: x
                        if isinstance(x, str)
                        else x * read_count_multiplier,
                    }
                }

        metric_match = re.search("^(?P<entity>bases|reads) in (?P<region>.+)", entity)
        if metric_match:
            entity = metric_match["entity"]
            REGION = metric_match["region"]
            region = improve_region(REGION)
            if entity == "bases":
                description = (
                    "Number ({}) of uniquely mapped bases to ".format(
                        config.base_count_desc
                    )
                    + region
                )
                configs = {
                    "min": 0,
                    "format": base_format,
                    "title": config.base_count_prefix + " Aln bases",
                    "description": description + ".",
                    "modify": lambda x: x
                    if isinstance(x, str)
                    else x * base_count_multiplier,
                }
                configs2 = {
                    "min": 0,
                    "max": 100,
                    "suffix": " %",
                    "format": base_format,
                    "title": "Aln bases on trg",
                    "description": description
                    + " relative to the number of uniquely mapped bases to the genome.",
                }
            else:
                description = (
                    "Number ({}) of uniquely mapped reads to ".format(
                        config.read_count_desc
                    )
                    + region
                    + ". DRAGEN V3.4 - V3.8:"
                    + " When region is the target BED, this metric is equivalent to and replaces"
                    + " Capture Specificity based on target region."
                    + " DRAGEN V3.9 - V4.1: Only reads with with MAPQ >= 1 are included."
                    + " Secondary and supplementary alignments are ignored."
                )
                configs = {
                    "min": 0,
                    "format": read_format,
                    "title": config.read_count_prefix + " Aln reads",
                    "description": description,
                    "modify": lambda x: x
                    if isinstance(x, str)
                    else x * read_count_multiplier,
                }
                configs2 = {
                    "min": 0,
                    "max": 100,
                    "suffix": " %",
                    "format": read_format,
                    "description": "Number of uniquely mapped reads to "
                    + region
                    + " relative to the number of uniquely mapped reads to the genome.",
                    "title": "Aln " + entity + " on trg",
                }

            return {"configs": configs, "configs2": configs2, "region": REGION}

        # Else unrecognized.
        metric = metric_pattern_match.string
        return {
            "configs": {
                "title": "Aln " + entity,
                "description": metric,
                "scale": "Purples",
            },
            "configs2": {
                "title": "Aln " + entity + " V2",
                "description": metric + ". Second value.",
                "scale": "Reds",
            },
            "warning": 1,
            "region": extract_region(metric),
        }

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
                "Number of uniquely mapped bases to "
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
                    "title": "Average " + entity + " V2",
                    "description": metric_pattern_match.string + ". Second value.",
                    "scale": "Reds",
                },
                "warning": 1,
                "region": REGION,
            }
        return {"configs": configs, "region": REGION}

    def get_Uniformity_configs(metric_pattern_match):
        """The following metrics are supported:
        Uniformity of coverage (PCT > float*mean) over region
        """
        entity = metric_pattern_match["entity"]
        REGION = metric_pattern_match["region"]
        region = improve_region(REGION)

        entity_match = re.search(r"\(pct\s*>\s*(\d+\.\d+)\*mean\)", entity)
        if entity_match:
            multiplier = entity_match.group(1)
            percent = str(float(multiplier) * 100) + "%"
            configs = {
                "suffix": " %",
                "min": 0,
                "max": 100,
                "title": ">" + multiplier + "*mean",
                "description": "Percentage of sites with coverage greater than "
                + percent
                + " of the mean coverage in "
                + region
                + ".",
            }
            if multiplier == "0.2":
                configs[
                    "description"
                ] += " Demonstrates the uniformity of coverage, the higher the better."

            return {"configs": configs, "region": REGION}

        # Else unrecognized.
        return {
            "configs": {
                "title": "Uniformity " + entity,
                "description": metric_pattern_match.string,
                "scale": "Purples",
            },
            "configs2": {
                "title": "Uniformity " + entity + " V2",
                "description": metric_pattern_match.string + ". Second value.",
                "scale": "Reds",
            },
            "warning": 1,
            "region": REGION,
        }

    # Special regex to match the [ix, inf) and [ix, jx)
    IX_JX = re.compile(r"\[(\d+)x:(inf|\d+)x?\)")

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
                title = GTQ + IX + "x"
                # Add extra entities to widen the title.
                if len(IX) < 6:
                    title += "&nbsp;" * (6 - len(IX))

                if IX == "0":
                    description += " with any coverage."
                else:
                    description += " with at least " + IX + "x coverage."
            else:
                title = entity
                if IX == "0" and JX == "1":
                    description += " with no coverage."
                else:
                    description += (
                        " with at least " + IX + "x but less than " + JX + "x coverage."
                    )
            return {
                "configs": {
                    "min": 0,
                    "max": 100,
                    "suffix": " %",
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
                "title": "PCT " + entity + " V2",
                "description": metric_pattern_match.string + ". Second value.",
                "scale": "Reds",
            },
            "warning": 1,
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
                "title": "Med " + entity + " V2",
                "description": metric_pattern_match.string + ". Second value.",
                "scale": "Reds",
            },
            "warning": 1,
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
                    "title": "Mean/Med aut cov",
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
                    "title": entity + " V2",
                    "description": metric_pattern_match.string + ". Second value.",
                    "scale": "Reds",
                },
                "warning": 1,
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
                "title": metric + " V2",
                "description": metric + " value2",
                "scale": "Reds",
            },
            "warning": 1,
            "region": extract_region(metric),
        }

    ALN_PAT["make_configs"] = get_Aligned_configs
    AVG_PAT["make_configs"] = get_Average_configs
    PCT_PAT["make_configs"] = get_PCT_configs
    UNI_PAT["make_configs"] = get_Uniformity_configs
    MED_PAT["make_configs"] = get_Median_configs
    RAT_PAT["make_configs"] = get_ratio_over_configs
    ANY_PAT["make_configs"] = get_any_configs

    # Finally return the closures.
    return coverage_metrics_parser, make_log_report


cov_parser, make_parsing_log_report = construct_coverage_parser()
