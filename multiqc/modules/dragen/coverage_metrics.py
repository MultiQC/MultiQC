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
from collections import defaultdict, OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file

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
    "title": None,           # Short title, table column title
    "description": None,     # Longer description, goes in mouse hover text
    "max": None,             # Minimum value in range, for bar / colour coding
    "min": 0,                # Maximum value in range, for bar / colour coding
    "ceiling": None,         # Maximum value for automatic bar limit
    "floor": None,           # Minimum value for automatic bar limit
    "minRange": None,        # Minimum range for automatic bar
    "scale": "GnBu",         # Colour scale for colour coding. False to disable.
    "bgcols": None,          # Dict with values: background colours for categorical data.
    "colour": "15, 150, 255",# Colour for column grouping
    "suffix": "",            # Suffix for value (eg. "%")
    "format": "{:,.2f}",     # Value format string - default 2 decimal places
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,# Styles for conditional formatting of table cell values
    "shared_key": None,      # See the link for description
    "modify": None,          # Lambda function to modify values, special case, see below
    "hidden": True,          # Set to True to hide the column in the general table on page load.
    "hidden_own": False,     # For non-general plots in own coverage sections.
    #"exclude": True         # Exclude all headers from the general html table.
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
    "hidden_own": False,   # For non-general plots in the coverage section.
    "exclude": False,      # True to exclude metric from the general html table.
    "exclude_own": False,  # True to exclude from own coverage section html table.
    "order_priority": None,# Used to specify columns' order in all tables.
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
BED = "target bed" # May be incorrect. No real data was available.

# The REGIONS holds the variables. You can add common region's settings here.
# These configs have higher priority than those in the SINGLE_HEADER.
# They also overwrite the automatically created ones.
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
  Otherwise long values can overlap with adjacent values in right columns.
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
  to none. "Configure Columns" is not included in general table's div block.

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

- Most configs below are for demonstration purposes only.
  You can adapt them according to your personal needs.

- You have the final word in defining settings.
  The code trusts you and does not check the chosen values.
"""

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
    # "aligned bases" has no support for region-specificity.
    # So, you can not set WGS/QC/BED-specific settings here.
    "aligned bases": {
        "title": "Aln bases",
        "scale": "RdYlGn",
        "colour": "0, 0, 255",
    },
    "aligned bases in region": {
        "title": "Bases on target",
        "scale": "Reds",
        "colour": "0, 0, 255",
        WGS: {
            "scale": "Purples",
        },
        QC: {
            "scale": "Greens",
        },
    },
    "aligned bases in region" + V2: {
        "title": "Bases on trg pct",
        "max": 100,
        "suffix": " %",
        "bgcols": {"NA": "#00FFFF",},
    },
    # "aligned reads" has no support for region-specificity.
    # So, you can not set WGS/QC/BED-specific settings here.
    "aligned reads": {
        "title": "Aln reads",
    },
    "aligned reads in region": {
        "title": "Reads on target",
        "scale": "RdGy",
        "colour": "255, 0, 0",
    },
    "aligned reads in region" + V2: {
        "title": "Reads on trg pct",
        "max": 100,
        "suffix": " %",
        "scale": "RdGy",
        "colour": "255, 0, 0",
    },
    "average alignment coverage over region": {
        "order_priority": 3,
        "title": "Depth",
        "scale": "BrBG",
        "colour": "0, 255, 255",
    },
    "average autosomal coverage over region": {
        "title": "Mean aut cov",
        "suffix": " x",
        "bgcols": {"NA": "#FF00FF"},
        WGS: {
            "description": """
            Average autosomal coverage over genome. Calculated as the
            number of bases that aligned to the autosomal loci in genome
            divided by the total number of loci in the autosomal loci in
            genome. If there are no autosomes in the reference genome, or
            the region does not intersect autosomes, this metric shows as NA.
            """,
        },
    },
    "average mitochondrial coverage over region": {
        "title": "MT cov&nbsp;&nbsp;",
        "scale": "Reds",
        "suffix": " x",
        "colour": "255, 0, 255",
        "order_priority": 4,
    },
    "average chr x coverage over region": {
        "title": "X cov",
        "suffix": " x",
        "colour": "255, 255, 0",
        "order_priority": 1,
    },
    "average chr y coverage over region": {
        "title": "Y cov",
        "suffix": " x",
        "colour": "255, 255, 0",
        "order_priority": 2,
    },
    # This metric can handle any float. "0.0" is just a technical trick.
    # Only the ">" operator and the "*mean" are supported.
    # Title and description shall not be set here, they are automated.
    # But if you want to handle a certain number manually, then you can
    # add an extra number-string key with desired configs.
    "uniformity of coverage (pct > 0.0*mean) over region": {
        "suffix": " %",
        "colour": "55, 255, 55",
        "hidden": False,
        "extra": {
            # "Uniformity of coverage (PCT > 0.2*mean) over region" metric
            "0.2": {
                "scale": "PiYG",
                "order_priority": 5,
                QC: {
                    "description": """Percentage of sites with coverage greater than
                    20% of the mean coverage in QC coverage region. Demonstrates
                    the uniformity of coverage, the higher the better.""",
                },
                # "Uniformity of coverage (PCT > 0.2*mean) over genome" metric
                WGS: {
                    "description": """Percentage of sites with coverage greater than
                    20% of the mean coverage in genome. Demonstrates
                    the uniformity of coverage, the higher the better.""",
                },
            },
            "0.4": {
                QC: {
                    "description": """Percentage of sites with coverage greater than
                    40% of the mean coverage in QC coverage region.""",
                },
            },
        },
    },
    # The following two metrics are the most tedious. The creation of the title
    # and description are automated. If you want to target some specific case,
    # then you can specify a certain combination of (ix, inf) or (ix, jx).
    # Just write a group in the "extra", you can consider it as a separate
    # metric, so you can also add general or region-specific configs there.
    # "0.x" are just a tecnical trick. Any natural number is meant instead.
    # "PCT of region with coverage [ix, inf)" metrics:
    "pct of region with coverage [0x:inf)": {
        "max": 100,
        "scale": "Blues",
        "suffix": " %",
        "colour": "255, 50, 25",
        WGS: {
            "scale": "Purples",
        },
        QC: {
            "scale": "Oranges",
        },
        "extra": {
            # "PCT of region with coverage [0x: inf)" metric.
            ("0x", "inf"): {
                "title": GTQ + "0x" + 8 * "&nbsp;",
                # Simple example of the special pair.
                "cond_formatting_rules": {
                    "red": [{"s_contains": ""}],# Each value is red/false by default
                    "green": [{"eq": 100.0}],   # If equal to 100, then green
                },
                "cond_formatting_colours": [
                    {"red": "#FF0000"},
                    {"green": "#00FF00"},
                ],
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
            ("1x", "inf"): {
                # "title": "Some title",
                # "description": "Some description", ...
            },
        },
    },
    # "PCT of region with coverage [ix, jx)" metrics:
    "pct of region with coverage [0x:0x)": {
        "max": 100,
        "scale": "Blues",
        "suffix": " %",
        "hidden_own": True,
        #"exclude": True,
        #"exclude_own": True,
        "extra": {
            ("0x", "1x"): {
                "title": "0x&nbsp;&nbsp;&nbsp;",
                WGS: {
                    "scale": "Oranges",
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
    "median autosomal coverage over region": {
        "title": "Med aut cov",
        "suffix": " x",
    },
    "mean/median autosomal coverage ratio over region": {
        "title": "Mean/med aut cov",
        "suffix": " x",
    },
    "xavgcov/yavgcov ratio over region": {
        "title": "XAvgCov/YAvgCov",
    },
    "xavgcov/autosomalavgcov ratio over region": {
        "title": "XAvgCov/AutAvgCov",
    },
    "yavgcov/autosomalavgcov ratio over region": {
        "title": "YAvgCov/AutAvgCov",
    },
}

# The TABLE_CONFIG defines common configs for all whole tables.
TABLE_CONFIG = {
    "namespace": NAMESPACE,# Name for grouping. Prepends desc and is in Config Columns modal
    "id": None,           # ID used for the table
    "table_title": None,  # Title of the table. Used in the column config modal
    "save_file": False,   # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sortRows": True,     # Whether to sort rows alphabetically
    "only_defined_headers": True, # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False, # Force a table to always be plotted (beeswarm by default if many rows)
}
# Below are region-specific configs with higher priority.
REGION_TABLE_CONFIG = {
    WGS: {
        "table_title": "WGS",
        # "no_beeswarm": True,
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
    """ Public members of the DragenCoverageMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined.
    Other methods can be added as well to provide extra features (eg module interface, JSON).
    """
    def add_coverage_metrics(self):
        """ The main function of the dragen coverage metrics module.
        The public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns either a view object containing
        a list of sample names, or an empty set if data could not be collected.
        """
        """ The general structure of the cov_data:

        cov_data
        |
        |_sample1___________________________________________________________________
        |         |         |                    |               |                  |
        |_sample2_|_________|____________________|_______________|__________________|
        |         |         |                    |               |                  |
        |_sampleN_|_________|____________________|_______________|__________________|
                 /          |                    |               |                  |
              __/           |                    |               |                  |
        genome     QC coverage region        target bed     target region    some future region
          |                 |                    |               |                  |
         wgs     qc-coverage-region-1,2...   target_bed     target_region    future-region-1,2...
          |                 |                    |               |                  |
        Metrics          Metrics              Metrics         Metrics         old/new Metrics
          |                 |                    |               |                  |
        Values           Values               Values          Values             Values

        First layer of each sample references the regions found in CSV files.
        The <coverage region prefix> are stored in the second layer.
        If region could not be extracted, then the "Unknown region" is used.
        """
        cov_data = defaultdict(lambda: defaultdict(dict))
        cov_headers = {}

        # Stores full path to properly match overall_mean.csv files later.
        orig_sample_names = defaultdict(lambda: defaultdict(dict))

        # Stores cleaned sample names with their filehandlers.
        # Used later to check for duplicates and inform if found.
        cleaned_sample_names = defaultdict(lambda: defaultdict(list))

        for file_handler in self.find_log_files("dragen/coverage_metrics"):
            out = cov_parser(file_handler, cov_headers)
            if out["success"]:
                self.add_data_source(file_handler, section = "stats")

                orig_s_name = out["sample_name"]
                clean_s_name = self.clean_s_name(orig_s_name, file_handler)
                phenotype = out["phenotype"]

                orig_sample_names[file_handler["root"]][orig_s_name][phenotype] = clean_s_name
                cleaned_sample_names[clean_s_name][phenotype].append(file_handler)

                # Add/overwrite the sample.
                cov_data[clean_s_name][out["region"]][phenotype] = out["data"]

        if not cov_data:
            return set()
        '''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
        Coverage data is collected now. The headers are also prepared.
        If other files (eg json) are needed, then these can be found
        with find_log_files, processed and used afterwards.
        See afterqc.py or supernova.py for example.
        CSS/JS can be imported, see fastqc.py for details.
        '''
        cov_data = self.ignore_samples(cov_data)

        # Check samples for duplicates.
        check_duplicate_samples(cleaned_sample_names)

        # Write data into the general table.
        gen_data, gen_headers = make_general_stats(cov_data, cov_headers)
        self.general_stats_addcols(gen_data, gen_headers, namespace = NAMESPACE)

        # Write general data to file.
        self.write_data_file(gen_data, "dragen_cov_metrics")

        # Extract coverage bed/target bed/wgs from _overall_mean_cov.csv files.
        # And prepare <coverage region prefix>-specific texts.
        bed_text = extract_source_files(self.overall_mean_cov_data, orig_sample_names)
        coverage_sections = make_cov_sections(cov_data, cov_headers, bed_text)

        # Special closure for reporting found info/warnings/errors,
        # which were collected while calling the cov_parser.
        # You can disable it anytime, if it is not wanted.
        make_parsing_log_report()

        for cov_section in coverage_sections:
            self.add_section(
                name = cov_section["name"],
                anchor = cov_section["anchor"],
                description = cov_section["description"],
                helptext = cov_section["helptext"],
                plot = cov_section["plot"],
            )
        return cov_data.keys()


def create_table_handlers():
    """Create data/headers for tables.
    Output: 2 closures
    1) make_general_stats
    2) make_own_coverage_sections
    """
    def order_headers(headers):
        """ Inserts single headers in the specified order."""
        indexes = set()
        ordered_headers = defaultdict(list)
        not_ordered_headers = {}
        for header in headers:
            if "order_priority" in headers[header]:
                order_priority = headers[header]["order_priority"]
                if isinstance(order_priority, int) or isinstance(order_priority, float):
                    indexes.add(order_priority)
                    ordered_headers[order_priority].append((header, headers[header]))
                else:
                    not_ordered_headers[header] = headers[header].copy()
            else:
                not_ordered_headers[header] = headers[header].copy()

        # Convert to list and sort in ascending order.
        if indexes:
            indexes = list(indexes)
            indexes.sort()
        else:
            return headers

        output_headers = OrderedDict()
        for index in indexes:
            for header in ordered_headers[index]:
                output_headers[header[0]] = {
                    config: val
                    for config, val in header[1].items()
                    if config not in EXTRA_HEADER
                }
        output_headers.update(not_ordered_headers)
        return output_headers

    def clean_region_prefix(phenotype, mode = ""):
        """Change <coverage region prefix> to make the output table better looking.
        Input:  phenotype, that is the <coverage region prefix> in the official standard.
        Output: a cleaned version of the phenotype.
        """
        phenotype = phenotype.strip()
        if re.search("QC", phenotype, re.IGNORECASE):
            phenotype = re.sub("qc", "QC", phenotype, flags = re.IGNORECASE)
            if mode == "Table header":
                phenotype = re.sub(r"-|_|\.", " ", phenotype)
            elif mode == "Column header":
                phenotype = re.sub(r"-|_|\.", "", phenotype)
                phenotype = re.sub("coverage|region", "", phenotype, flags = re.IGNORECASE)
            return phenotype
        if re.search("wgs", phenotype, re.IGNORECASE):
            phenotype = phenotype.upper()

        phenotype = re.sub(r"-|_|\.", " ", phenotype)
        return phenotype

    def set_table_config(region):
        table_config = {config: val for config, val in TABLE_CONFIG.items() if val is not None}
        if region in REGION_TABLE_CONFIG and REGION_TABLE_CONFIG[region]:
            table_config.update(REGION_TABLE_CONFIG[region])
        return table_config

    def make_general_stats(data, headers):
        """Prepare data and headers for the general stats table."""
        gen_data = defaultdict(dict)
        gen_headers = defaultdict(dict)
        for sample in data:
            for region in data[sample]:
                for phenotype in data[sample][region]:
                    pheno = data[sample][region][phenotype]
                    for metric in pheno:
                        if "exclude" in headers[metric] and headers[metric]["exclude"]:
                            pass
                        else:
                            new_metric_id = metric + "-" + phenotype
                            gen_data[sample][new_metric_id] = pheno[metric]
                            gen_headers[new_metric_id] = headers[metric].copy()
                            subregion_cleaned = clean_region_prefix(phenotype, mode = "Column header")
                            if "title" in gen_headers[new_metric_id]:
                                gen_headers[new_metric_id]["title"] = (
                                    re.sub("&nbsp;", "", gen_headers[new_metric_id]["title"])
                                    + " " + subregion_cleaned
                                )
        return gen_data, order_headers(gen_headers)


    def make_own_coverage_sections(data_input, headers_input, bed_textes):
        """Create non-general region-specific sections."""
        plots = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
        for sample in data_input:
            for region in data_input[sample]:
                reg = data_input[sample][region]
                for phenotype in reg:
                    data = {}
                    headers = defaultdict(dict)
                    for metric in reg[phenotype]:
                        if "exclude_own" in headers_input[metric] and \
                            headers_input[metric]["exclude_own"]:
                            pass
                        else:
                            id_metric = metric + "-own-table-" + phenotype
                            data[id_metric] = reg[phenotype][metric]
                            headers[id_metric] = headers_input[metric].copy()

                            if "hidden_own" in headers[id_metric]:
                                headers[id_metric]["hidden"] = headers[id_metric]["hidden_own"]

                    plots[region][phenotype]["data"][sample] = data
                    plots[region][phenotype]["headers"] = order_headers(headers)
                    plots[region][phenotype]["config"] = set_table_config(region)

        sections = []
        for region in plots:
            for phenotype in plots[region]:
                cleaned_region = clean_region_prefix(phenotype, mode = "Table header")
                sections.append({
                    "name": cleaned_region + " Metrics",
                    # Spaces are replaced with hyphens. MultiQC linter complains otherwise.
                    "anchor": "dragen-cov-metrics-own-sec-" + re.sub("\s+", "-", region) + re.sub("\s+", "-", phenotype),
                    "helptext": """
                    The following criteria are used when calculating coverage:\n
                    * Duplicate reads and clipped bases are ignored.\n
                    * Only reads with `MAPQ` > `min MAPQ` and bases with `BQ` > `min BQ` are considered\n
                    Considering only bases usable for variant calling, _i.e._ excluding:\n
                    1. Clipped bases\n
                    2. Bases in duplicate reads\n
                    3. Reads with `MAPQ` < `min MAPQ` (default `20`)\n
                    4. Bases with `BQ` < `min BQ` (default `10`)\n
                    5. Reads with `MAPQ` = `0` (multimappers)\n
                    6. Overlapping mates are double-counted\n
                    """,
                    "description": """
                    Coverage metrics over a {reg}. {source_file}. Press the `Help` button for details.\n
                    """.format(reg = region, source_file = bed_textes[phenotype]),
                    "plot": table.plot(plots[region][phenotype]["data"],
                                       plots[region][phenotype]["headers"],
                                       plots[region][phenotype]["config"])
                })
        return sections

    # Return the closures.
    return make_general_stats, make_own_coverage_sections
make_general_stats, make_cov_sections = create_table_handlers()


def check_duplicate_samples(sample_names):
    """ Check samples for duplicate names. Warn about found ones."""
    message = ""
    for sample_name in sample_names:
        for phenotype in sample_names[sample_name]:
            if len(sample_names[sample_name][phenotype]) > 1:
                message += "  " + sample_name + "\n  was built from following samples:\n"
                for file_handler in sample_names[sample_name][phenotype]:
                    message += "    " + file_handler["root"] + ": " + file_handler["fn"] + "\n"
    if message:
        message = (
            "\n\nDuplicate sample names were found. " +
            "The last one overwrites previous data.\n\n" +
            message + "\n"
        )
    log.warning(message)


def extract_source_files(overall_mean_cov_data, sample_names):
    """Matches _overall_mean_cov.csv to corresponding _coverage_metrics.csv
    Extracts bed file names and creates
    Returns a dict with phenotype: [bed source files].
    The values can be also checked.
    """
    phenotypes_with_data = defaultdict(lambda: defaultdict(list))

    for root in overall_mean_cov_data:
        for orig_sample in overall_mean_cov_data[root]:
            for phenotype in overall_mean_cov_data[root][orig_sample]:
                if (root in sample_names) and \
                   (orig_sample in sample_names[root]) and \
                   (phenotype in sample_names[root][orig_sample]):
                    data = overall_mean_cov_data[root][orig_sample][phenotype]
                    clean_sample = sample_names[root][orig_sample][phenotype]
                    phenotypes_with_data[phenotype][data["source_file"]].append(clean_sample)

    beds_for_sections = {}
    for phenotype in phenotypes_with_data:
        text = ""
        bed_sources = list(phenotypes_with_data[phenotype].keys())
        # The lenght shall be one, but who knows.
        if len(bed_sources) == 1:
            text = "All samples are based on the " + bed_sources[0].split("/")[-1]
        # If not, then iterate through source files.
        else:
            for source in bed_sources:
                text += "The " + source + " is based on:"
                for sample in phenotypes_with_data[phenotype][source]:
                    text += sample + ", "

        beds_for_sections[phenotype] = text

    return beds_for_sections


def construct_coverage_parser():
    """ Isolation for all parsing machinery.
    Returns 2 closures:
    1) Parser for coverage data stored in csv files.
    2) Log reporter for found info/warnings/errors.
    """
    # All regexes are constructed to be as general and simple as possible to speed up the matching.
    # "make_configs" will point later to a function, which automatically sets some header's configs.
    ALN_PAT = {
        "RGX": re.compile("^Aligned (?P<entity>.+)", re.IGNORECASE),
        "make_configs": None
    }
    AVG_PAT = {
        "RGX": re.compile("^Average (?P<entity>.+?) coverage over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None
    }
    PCT_PAT = {
        "RGX": re.compile("^PCT of (?P<region>.+?) with coverage (?P<entity>.+)", re.IGNORECASE),
        "make_configs": None
    }
    UNI_PAT = {
        "RGX": re.compile("^Uniformity of coverage (?P<entity>.+?) over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None
    }
    MED_PAT = {
        "RGX": re.compile("^Median (?P<entity>.+?) coverage over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None
    }
    RAT_PAT = {
        "RGX": re.compile("(?P<entity>.+?) ratio over (?P<region>.+)", re.IGNORECASE),
        "make_configs": None
    }
    ANY_PAT = {
        "RGX": re.compile(".+", re.IGNORECASE),
        "make_configs": None
    }
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
    If some new metric-specific peculiarities arise, they can be added to the make_metric_id.
    '''
    def make_metric_id(metric):
        """ Tries to fix consistency issues that might arise in coverage metrics data."""

        # Single backslashes are escaped.
        metric_id = re.sub("\s+", " ", metric).strip().lower()

        pct_case = PCT_PAT["RGX"].search(metric_id)
        if pct_case:
            metric_id = "pct of " + pct_case["region"] + " with coverage " + pct_case["entity"].replace(" ", "")

        return metric_id


    def get_std_configs(configs_dict):
        """ Copies the standard real/virtual configurations."""
        return {
            config: value
            for config, value in configs_dict.items()
            if (config in SINGLE_HEADER or config in EXTRA_HEADER) and value is not None
        }

    FILE_RGX = re.compile(r"^([^.]+)\.(.+?)_coverage_metrics(.*).csv")
    LINE_RGX = re.compile("^COVERAGE SUMMARY,,([^,]+),([^,]+),?([^,]*)")

    # Used to simplify the code and speed up execution.
    cov_headers_support = defaultdict(dict)

    def coverage_metrics_parser(file_handler, headers):
        """ Parser for coverage metrics csv files.
        Input:  file_handler with necessary info - file name/content/root
        Output: {"success": 0} if file name test failed. Otherwise:
                {"success": 0/1, "sample_name": <output prefix>,
                 "phenotype": <coverage region prefix>,
                 "region": extracted_region,
                 "data": extracted_metrics_and_values
                } 0 if all file's lines are invalid.
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
                metric, value1, value2 = line_match.group(1), line_match.group(2), line_match.group(3)
                success = 1

            # Otherwise check if line is empty. If not then report it and go to the next line.
            else:
                if not re.search("^\s*$", line):
                    log_data["invalid_file_lines"][root_name][file_name].append(line)
                continue

            """ Check the extracted values.
            These are int/float and in some cases NA for now. All other values will be reported.
            The type conversion is only performed to int and float. Other types of data are left
            unchanged as strings. It has some consequences for the "modify" configuration.
            Details are in the description for the "make_auto_configs". No additional try blocks
            are required to support other simple types of values (eg Inf, None). Just adapt the regexes.
            """
            try:
                value1 = int(value1)
            except ValueError:
                try:
                    value1 = float(value1)
                except ValueError:
                    if not re.search("^NA$", value1, re.IGNORECASE):
                        log_data["unusual_values"][root_name][file_name][metric + V2] = value1

            if value2:
                try:
                    value2 = float(value2)
                except ValueError:
                    try:
                        value2 = int(value2)
                    except ValueError:
                        if not re.search("^NA$", value2, re.IGNORECASE):
                            log_data["unusual_values"][root_name][file_name][metric + V2] = value2

            metric_id = make_metric_id(metric)
            if metric_id not in headers:
                for MPAT in METRIC_PATTERNS:
                    metric_match = MPAT["RGX"].search(metric_id)
                    if metric_match:
                        output = MPAT["make_configs"](metric_match)

                        if "region" in output:
                            region = output["region"]
                            cov_headers_support[metric_id]["region"] = region

                        if "warning" in output:
                            cov_headers_support[metric_id]["warning"] = 1

                        configs = get_std_configs(SINGLE_HEADER)

                        if AUTO_CONFIGS_ENABLED:
                            auto_configs = output["configs"]
                            configs.update(auto_configs)

                        if USER_CONFIGS_ENABLED:
                            user_configs = make_user_configs(metric_id, region)
                            configs.update(user_configs)

                        headers[metric_id] = configs

                        if value2:
                            configs = get_std_configs(SINGLE_HEADER)

                            if AUTO_CONFIGS_ENABLED:
                                if "configs2" in output:
                                    auto_configs = output["configs2"]
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

            # In case a second value found, but the header was not filled
            if value2 and (metric_id + V2 not in headers):
                headers[metric_id + V2] = {
                    "title": metric,
                    "description": "Unrecognized metric.",
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
            "region": region,
            "phenotype": phenotype,
            "data": data,
        }

    def make_log_report():
        """The only purpose of this function is to create a readable and informative log output."""

        if log_data["invalid_file_names"]:
            log_message = ("\n\nThe file names must conform to the following structure:\n" +
                           "<output prefix>.<coverage region prefix>_coverage_metrics<arbitrary suffix>.csv\n\n" +
                           "The following files are not valid:\n"
            )
            for root in log_data["invalid_file_names"]:
                log_message += "  " + root + ":\n"
                for file in log_data["invalid_file_names"][root]:
                    log_message += "    " + file + "\n"

            log.warning(log_message + "\n")

        if log_data["invalid_file_lines"]:
            log_message = ("\n\nThe lines in files must be:\n" +
                           "COVERAGE SUMMARY,,<metric>,<value1> or " +
                           "COVERAGE SUMMARY,,<metric>,<value1>,<value2>\n\n"  +
                           "The following files contain invalid lines:\n"
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
            log_message = ("\n\nAll metrics' values except int, float and NA are non-standard.\n"  +
                           "The following files contain non-standard values:\n"
            )
            for root in log_data["unusual_values"]:
                log_message += "  " + root + ":\n"
                for file in log_data["unusual_values"][root]:
                    log_message += "    " + file + ":\n"
                    for metric in log_data["unusual_values"][root][file]:
                        log_message += ("      " + metric + " = " +
                                        log_data["unusual_values"][root][file][metric] + "\n")

            log.warning(log_message + "\n")


    def make_user_configs(metric, region):
        """Creates the user-defined configurations."""
        configs = {}

        # Check if empty, in REGIONS and defined.
        if region and region in REGIONS and REGIONS[region]:
            metric = metric.replace(region, "region")
            configs.update(get_std_configs(REGIONS[region]))

        pct_case = PCT_PAT["RGX"].search(metric)
        uni_case = UNI_PAT["RGX"].search(metric)
        if pct_case or uni_case:
            metric = re.sub("\d+", "0", metric)

        # Now set region-specific parameters for a given metric.
        if metric in METRICS:
            metric = METRICS[metric]
            configs.update(get_std_configs(metric))

            # If region exists and not empty.
            if region and region in metric:
                configs.update(get_std_configs(metric[region]))

        # Try to set (ix:inf)/(ix:jx)-specific settings.
        if pct_case and "extra" in metric and metric["extra"]:

            ix_jx = re.search(r"\[(\d+x):(inf|\d+x)\)", pct_case["entity"])
            if ix_jx:
                _group = (ix_jx.group(1), ix_jx.group(2))

                if _group in metric["extra"]:
                    metric = metric["extra"][_group]
                    configs.update(get_std_configs(metric))
                    if region in metric:
                        configs.update(get_std_configs(metric[region]))

        # Try to match the float.
        if uni_case and "extra" in metric and metric["extra"]:
            entity_match = re.search(r"\(PCT > (\d+\.\d+)\*mean\)", uni_case["entity"], re.IGNORECASE)
            if entity_match:
                _float = entity_match.group(1)
                if _float in metric["extra"]:
                    metric = metric["extra"][_float]
                    configs.update(get_std_configs(metric))
                    if region in metric:
                        configs.update(get_std_configs(metric[region]))
        return configs


    '''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    The following functions are mainly used to create configurations for headers.

    Input:
    - match_object, which stores the result of matching a metric.

    Output:
    A dictionary with keys:
    - obligatory  configs   dict with real/logical configs from the SINGLE_HEADER.
    - obligatory  configs2  if the second value is presented or
                            in case when metric can not be recognized.
    - optional    region    if presented in a metric and can be extracted easily.
                            Returning it is encouraged.
    - optional    warning   can be returned if metric could not be recognized.
                            The associated value is irrelevant.

    Please note: modify lambda/func must check for string input.
    """""""""""""""""""""""""""""""""""""""""""""""""""""""""'''

    ##---------------------------------------------------##
    ## This little segment was taken from ./dragen/utils.py
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
    ##---------------------------------------------------##

    """
    def make_hex(red, green, blue):
        "Returns a color int to hex str."
        red = min(abs(red), 255)
        red = ('0' if red < 16 else "") + hex(red)[2:]
        green = min(abs(green), 255)
        green = ('0' if green < 16 else "") + hex(green)[2:]
        blue = min(abs(blue), 255)
        blue = ('0' if blue < 16 else "") + hex(blue)[2:]
        return "#" + red + green + blue

    Different distributions can be used: uniform, exp, Poisson, Gauss, gamma, Rayleigh...
    It is not even necessary to implement them, just create a chunk of data in advance
    and store it as variable here.

    Red_to_Green_unif_on100 = [
        {str(color_n): make_hex(color_n, 255-color_n, 0)}
        for color_n in range(256)
    ].append({"1": "#FF00FF"})

    Red_to_Green_unif_on100 = {
        {str(color_n): make_hex(color_n, 255-color_n, 0)}
        for color_n in range(256)
    }.update({"magenta": "#FF00FF"})
    f_rules = {"red": [{"s_contains": ""}], "magenta": [{"gt": 100}]}
    f_rules = {[{"gt": ((green)/255)*100}]}

    {"red": [{"s_contains": ""}], "magenta": [{"gt": 100}]}
    for color in range(256):
        colour = "color" + str(green)
        f_rules[colour] =
    configs["cond_formatting_rules"] = f_rules
    configs["cond_formatting_colours"] = f_colours
    """

    def get_Aligned_configs(match_object):
        """The following metrics are supported:
        Aligned bases, Aligned reads
        Aligned bases in region, Aligned reads in region
        """
        configs = {}
        entity = match_object["entity"].lower()

        metric_match = re.search("^(bases|reads)$", entity)
        if metric_match:
            count_desc = config.base_count_desc if entity == "bases" else config.read_count_desc
            description = "Total number ({}) of aligned ".format(count_desc) + entity + "."
            title = (config.base_count_prefix if entity == "bases" else config.read_count_prefix) + " Aln " + entity
            base_or_read_format = base_format if entity == "bases" else read_format
            count_multiplier = base_count_multiplier if entity == "bases" else read_count_multiplier
            modify = lambda x: x if isinstance(x, str) else x * count_multiplier
            configs.update({
                "min": 0, "format": base_or_read_format,
                "description": description,
                "title": title,
                "modify": modify,
            })
            return {"configs": configs}

        configs2 = {}
        metric_match = re.search("^(?P<entity>bases|reads) in (?P<region>.+)", match_object["entity"])
        if metric_match:
            entity = metric_match["entity"].lower()
            region = metric_match["region"]
            configs.update({
                "min": 0, "format": base_format if entity == "bases" else read_format,
                "title": (config.base_count_prefix if entity == "bases" else config.read_count_prefix) + " Aln " + entity,
                "description":
                    ("Number ({}) of uniquely mapped bases to ".format(config.base_count_desc) + region + ".") if entity == "bases" else
                     "Number ({}) of uniquely mapped reads to ".format(config.read_count_desc) + region +
                    ". When " + region + " is the target BED, this metric is equivalent to and replaces Capture Specificity based on target region."
            })
            if entity == "reads": configs["modify"] = lambda x: x if isinstance(x, str) else x * read_count_multiplier
            else:                 configs["modify"] = lambda x: x if isinstance(x, str) else x * base_count_multiplier

            configs2.update({"min": 0, "max": 100, "suffix": " %",
                            "format": base_format if entity == "bases" else read_format,
                            "description": "Percentage of aligned " + entity + " in " + region + ".",
                            "title": "Aln " + entity + " on trg"})

            return {"configs": configs, "configs2": configs2, "region": metric_match["region"]}

        # Else unrecognized.
        configs["title"] = "Aln " + entity
        configs2["title"] = "Aln " + entity + " value2"
        return {"configs": configs, "configs2": configs2, "warning": 1}

    def get_Average_configs(match_object):
        """The following metrics are supported:
        Average alignment coverage over region
        Average chromosome X coverage over region, Average chromosome Y coverage over region
        Average mitochondrial coverage over region
        Average autosomal coverage over region
        """
        configs = {}
        entity = match_object["entity"]
        region = match_object["region"]

        entity_match = re.search("^(alignment|chr x|chr y|mitochondrial|autosomal)$", entity)
        if entity_match:
            configs.update({"suffix": " x", "min": 0, "title": "Depth"})

            if entity == "alignment":
                configs["description"] = """
                Number of uniquely mapped bases to {region}
                divided by the number of sites in {region}.
                """.format(region = region)
                configs["title"] = "AvgAlnCov"

            elif entity == "autosomal":
                configs["description"] = """
                Total number of bases that aligned to the autosomal loci in {region}
                divided by the total number of loci in the autosomal loci in {region}.
                If there is no autosome in the reference genome, or the {region}
                does not intersect autosomes, this metric shows as NA.
                """.format(region = region)
                configs["title"] = "AvgAutCov"

            elif entity == "mitochondrial":
                configs["description"] = """
                Total number of bases that aligned to the intersection of the mitochondrial
                chromosome with {region} divided by the total number of loci in the
                intersection of the mitochondrial chromosome with {region}. If there is no
                mitochondrial chromosome in the reference genome or the {region} does not
                intersect mitochondrial chromosome, this metric shows as NA.
                """.format(region = region)
                configs["title"] = "AvgMitoCov"

            else:
                entity = re.search("(x|y)", entity).group(1).capitalize()
                configs["title"] = "Avg" + entity + "cov"
                entity = "chromosome " + entity
                configs["description"] = """
                Total number of bases that aligned to the intersection of {entity}
                with {region} divided by the total number of loci in the intersection of
                {entity} with {region}. If there is no {entity} in the reference genome or the
                {region} does not intersect {entity} this metric shows as NA.
                """.format(region = region, entity = entity)
            return {"configs": configs, "region": match_object["region"]}

        # Else unrecognized.
        entity = match_object["entity"]
        configs["title"] = "Avg " + entity
        configs2 = {"title": "Avg " + entity + " val2"}
        return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

    def get_Uniformity_configs(match_object):
        """The following metrics are supported:
        Uniformity of coverage (PCT > float*mean) over region
        """
        configs = {}

        entity = match_object["entity"]
        region = match_object["region"]


        entity_match = re.search(r"\(PCT\s*>\s*(\d+\.\d+)\*mean\)", entity, re.IGNORECASE)
        if entity_match:
            multiplier = entity_match.group(1)
            percent = str(float(multiplier)*100) + "%"

            configs.update({"suffix": " %", "min": 0, "max": 100, "title": ">" + multiplier + "*mean",
                            "description": "Percentage of sites with coverage greater than " +
                            percent + " of the mean coverage in " + region + ".",})
            if multiplier == "0.2":
                configs["description"] = """
                Percentage of sites with coverage greater than 20%
                of the mean coverage in {region}. Demonstrates the
                uniformity of coverage, the higher the better.
                """.format(region = region)

            return {"configs": configs, "region": match_object["region"]}

        # Else unrecognized.
        configs2 = {}
        configs["title"] = "Uniformity " + entity
        configs2["title"] = "Uniformity " + entity + " value2"
        return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

    def get_PCT_configs(match_object):
        """The following metrics are supported:
        PCT of region with coverage [ix, inf)
        PCT of region with coverage [ix, jx)
        """
        configs = {}
        entity = match_object["entity"]
        region = match_object["region"]

        entity_match = re.search(r"\[(\d+)x:(inf|\d+)x?\)", entity)
        if entity_match:
            left_ix  = entity_match.group(1)
            right_jx = entity_match.group(2)
            nbsp_entities = "&nbsp;"*(6 - len(left_ix))
            title = GTQ + left_ix + "x" + nbsp_entities if "inf" == right_jx else entity
            configs.update({
                "suffix": " %", "min": 0, "max": 100, "title": title,
                "description":(
                    "Percentage of sites in " + region + " with no coverage." if (left_ix == "0" and right_jx == "1") else
                   ("Percentage of sites in " + region + " with at least " + left_ix + "x coverage.")
                    if re.search("inf", right_jx, re.IGNORECASE) else
                   ("Percentage of sites in " + region + " with at least " + left_ix + "x but less than " + right_jx + "x coverage.")),
            })
            return {"configs": configs, "region": match_object["region"]}

        # Else unrecognized.
        configs2 = {}
        configs["title"] ="PCT " + entity
        configs2["title"] = "PCT " + entity + " value2"
        return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

    def get_Median_configs(match_object):
        """The following metrics are supported:
        Median autosomal coverage over region
        """
        configs = {}
        entity = match_object["entity"]
        region = match_object["region"]

        metric_match = re.search("^(autosomal)$", entity)
        if metric_match:
            configs.update({"suffix": " x", "min": 0, "title": "Med aut cov",
                            "description": "Median alignment coverage over the autosomal loci in " + region +
                                            ". If there is no autosome in the reference genome or the " + region +
                                            " does not intersect autosomes, this metric shows as NA.",
            })
            return {"configs": configs, "region": match_object["region"]}

        # Else unrecognized.
        configs2 = {}
        configs["title"] = "Med " + entity
        configs2["title"] = "Med " + entity + " value2"
        return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

    def get_ratio_configs(match_object):
        """The following metrics are supported:
        Mean/Median autosomal coverage ratio over region
        XAvgCov/YAvgCov ratio over region
        XAvgCov/AutosomalAvgCov ratio over region
        YAvgCov/AutosomalAvgCov ratio over region
        """
        configs = {}
        entity = match_object["entity"]
        region = match_object["region"]

        entity_match = re.search("^Mean/Median autosomal coverage$", entity, re.IGNORECASE)
        if entity_match:
            configs.update({"suffix": " x", "title": "Mean/med aut cov", "min": 0,
                            "description": "Mean autosomal coverage in " + region +
                            " divided by the median autosomal coverage in " + region +
                            ". If there is no autosome in the reference genome or the " + region +
                            " does not intersect autosomes, this metric shows as NA."
            })
            return {"configs": configs, "region": match_object["region"]}

        entity_match = re.search("^(XAvgCov|YAvgCov|AutosomalAvgCov)/(XAvgCov|YAvgCov|AutosomalAvgCov)$", entity, re.IGNORECASE)
        if entity_match:
            devident = re.search("(X|Y|Aut)", entity_match.group(1), re.IGNORECASE).group(1).capitalize()
            devisor  = re.search("(X|Y|Aut)", entity_match.group(2), re.IGNORECASE).group(1).capitalize()
            configs.update({"min": 0, "title": devident + "AvgCov/" + devisor + "AvgCov"})

            chrom = ("chromosome " + devident) if (len(devident) == 1) else "autosomal chromosome"
            chrom2  = ("chromosome " + devisor)  if (len(devisor)  == 1) else "autosomal chromosome"
            configs["description"] = ("Average " + chrom + " alignment coverage in " + region +
                                      " divided by the average " + chrom2 + " alignment coverage in  " + region)

            return {"configs": configs, "region": match_object["region"]}

        # Else unrecognized.
        configs2 = {}
        configs["title"] = entity
        configs2["title"] = entity + " value2"
        return {"configs": configs, "configs2": configs2, "warning": 1, "region": match_object["region"]}

    def get_any_configs(match_object):
        """This function handles all other metrics, whose general structure could not be recognized."""
        configs = {"title": match_object.string, "scale": "Purples",}
        configs2 = {"title":  match_object.string + " value2", "scale": "Reds",}
        return {"configs": configs, "configs2": configs2, "warning": 1}

    ALN_PAT["make_configs"] = get_Aligned_configs
    AVG_PAT["make_configs"] = get_Average_configs
    PCT_PAT["make_configs"] = get_PCT_configs
    UNI_PAT["make_configs"] = get_Uniformity_configs
    MED_PAT["make_configs"] = get_Median_configs
    RAT_PAT["make_configs"] = get_ratio_configs
    ANY_PAT["make_configs"] = get_any_configs

    # Finally return the closures.
    return coverage_metrics_parser, make_log_report
cov_parser, make_parsing_log_report = construct_coverage_parser()
