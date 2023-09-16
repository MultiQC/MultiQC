'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers ploidy estimation data and prepares it for the output report.
It relies on the following official sources:
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/PloidyEstimator_fDG.htm
Section "Ploidy Calling" p.177:
https://emea.support.illumina.com/content/dam/illumina-support/documents/documentation/software_documentation/dragen-bio-it/Illumina-DRAGEN-Bio-IT-Platform-User-Guide-1000000141465-00.pdf


You can define a closed interval [lower_limit, upper_limit] of acceptable values for the metrics:
1 median / Autosomal median, 2 median / Autosomal median, ...,  22 median / Autosomal median

Example: add multiqc_config.yaml file in your current working directory with this content:
ploidy_estimation_user_configs:
    lower_limit: 0.75
    upper_limit: 1.25
'''

import logging
import re
from collections import OrderedDict, defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

from .utils import check_duplicate_samples, flatten_4D_data, make_gen_stats_headers, make_log_report, make_own_headers

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN ploidy estimation"


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The SINGLE_HEADER is used to define common settings for all headers:
https://multiqc.info/docs/development/plots/#creating-a-table
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
    "colour": "0, 0, 255",  # Colour for column grouping
    "suffix": None,  # Suffix for value (eg. "%")
    "format": "{:,.2f}",  # Value format string - default 2 decimal places
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values
    "hidden": True,  # Set to True to hide the column in the general table on page load.
    "hidden_own": True,  # For non-general plots in the ploidy section.
    "exclude": True,  # True to exclude all headers from the general html table.
}
'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER. These are just extensions,
but can be used as if they were a part of the SINGLE_HEADER. Extra configurations are not
set by default. Only the SINGLE_HEADER is used as a basis for a header. But you can always
append an extra config(eg hidden_own, exclude) to the SINGLE_HEADER to overcome this issue.

Hints and rules:
- "hidden" is used for own ploidy section if "hidden_own" is not provided.

- "exclude" and "exclude_own" may be useful in cases where html size matters.
  Eg. many samples but only few metrics of interest.

- "order_priority" can be used to arrange table's columns. Only int and float are valid.
'''
EXTRA_HEADER = {
    "hidden_own": None,  # For non-general plots in the ploidy section.
    "exclude": None,  # True to exclude metric from the general html table.
    "exclude_own": None,  # True to exclude from own ploidy section html table.
    "order_priority": None,  # Used to specify columns' order in all tables.
}

'''"""""""""""""""""""""""""""""""""""""""""""""""
Officially defined file format could not be found.
The links above show only the general non-technical example.
According to examined real data, there are 4 columns in total:
Section    Void-column    Metric    Value

The only section found in data was 'PLOIDY ESTIMATION'
It is unknown if it is the only valid section or others already exist/will be added.
The second column does not contain anything.
It is unknown if it is the final constant structure or it was reserved for future use.

Due to this uncertainty the following structure was chosen to store parsed ploidy data:
{sample: {section: {Region_or_Sample: {metric: value}}}}

Might look overcomplicated, but it conforms to the standard.

If you wish to add new sections, take a look at the make_consistent_section.
'''

## Sections:
PLOIDY_ESTIMATION = "PLOIDY ESTIMATION"
ADDITIONAL_SECTION = "NON-STANDARD ADDITIONAL SECTION"  # For EXTRA_METRICS.

## RG/Sample aka Second column:
EMPTY = ""
WGS = "WGS"  # Not detected/used. Just an example.

## Special case. Token for any section/region/sample.
ANY = object()

"""
Only INCLUDED_SECTIONS will be shown in html table.
The parser supports only the first value (4th column).
Can be used as a quick but unreliable way to control html table.
"""
INCLUDED_SECTIONS = {
    PLOIDY_ESTIMATION: [EMPTY],
    ADDITIONAL_SECTION: [ANY],
}

# INCLUDED_SECTIONS = {ANY:[ANY]}  # Include all sections/regions/samples.

# Used to detect and report what is not properly supported by the module.
# Works as INCLUDED_SECTIONS.
SUPPORTED_SECTIONS = {
    PLOIDY_ESTIMATION: [EMPTY],
}

'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The STD_METRICS is the container for abstract representations of the standard metrics.
They are designed to be seen as simple identifiers of the well known metrics.
You can set any real/virtual configuration defined in the SINGLE_HEADER.

Some rules:

- "modify" must conform to the following general structure:
  lambda x: (arbitrary expression) if isinstance(x, str) else (math expression)
  The x is either number or string. The latter is some value, which could not be
  converted to int or float. NA for instance. Such values can be modified in the
  (arbitrary expression) if needed. If not then just return the input.

- If you wish to add new metrics, take a look at the make_consistent_metric.
'''

RESERVED = 100  # Used for automatic ordering. Do not exceed it.

STD_METRICS = {
    "ploidy estimation": {
        "order_priority": 1,
        "exclude": False,
        "hidden": False,
        "hidden_own": False,
        "title": "Sex",
        "description": "Sex chromosome ploidy estimation (XX, XY, X0, 00, etc.)",
        "scale": "Set3",
        "colour": "255, 0, 0",
        "cond_formatting_rules": {
            # Each value is red/false by default
            "red": [{"s_contains": ""}],
            "green": [{"s_eq": "XX"}, {"s_eq": "XY"}, {"s_eq": "YX"}],
        },
        "cond_formatting_colours": [
            {"red": "#FF0000"},
            {"green": "#00FF00"},
        ],
    },
    "autosomal median coverage": {
        "order_priority": 3,
        "hidden_own": False,
        "title": "AutMedCov",
        "scale": "YlOrRd",
        "colour": "0, 0, 255",
        "description": "Autosomal median coverage.",
    },
    "x median coverage": {
        "order_priority": 4,
        "hidden_own": False,
        "title": "XMedCov",
        "scale": "BrBG",
        "colour": "0, 0, 255",
        "description": "X median coverage.",
    },
    "y median coverage": {
        "order_priority": 5,
        "hidden_own": False,
        "title": "YMedCov",
        "scale": "Purples",
        "colour": "0, 0, 255",
        "description": "Y median coverage.",
    },
    "x median / autosomal median": {
        "order_priority": 6,
        "title": "XMed/AutMed",
        "scale": "Greens",
        "colour": "0, 0, 255",
        "description": "X median coverage devided by autosomal median coverage.",
    },
    "y median / autosomal median": {
        "order_priority": 7,
        "title": "YMed/AutMed",
        "scale": "Oranges",
        "colour": "0, 0, 255",
        "description": "Y median coverage devided by autosomal median coverage.",
    },
    # The set of "1,2,3,...,22 median / Autosomal median" metrics.
    # title, description and order_priority are created automatically.
    # Each single metric can be redefined in the "extra".
    "autosome_n median / autosomal median": {
        "hidden_own": True,
        # "exclude_own": True,
        "colour": "255, 255, 0",
        "suffix": " %",
        "extra": {
            # "1 median / Autosomal median" metric:
            "1": {
                "description": "Median coverage of autosome 1 devided by autosomal median coverage.",
            },
            # "22 median / Autosomal median" metric:
            "22": {
                "description": "Median coverage of autosome 22 devided by autosomal median coverage.",
            },
        },
    },
}


# User/Default settings for the EXTRA_METRICS.
PLOIDY_MOD = getattr(config, "ploidy_estimation_user_configs", {})
LOWER_LIMIT = PLOIDY_MOD.get("lower_limit", 0.95)
UPPER_LIMIT = PLOIDY_MOD.get("upper_limit", 1.05)

# The EXTRA_METRICS contains non-standard "metrics".
# It is basically just a holder for some calculated results.
EXTRA_METRICS = {
    # Autosomes from "1,2,...,22 median / Autosomal median" metrics,
    # whose values are out of user-defined/default valid range.
    "autosome median to autosomal median ratios out of range": {
        "order_priority": 2,
        "hidden_own": False,
        # "title" is wide to avoid overlapping of adjacent rows.
        "title": "Autosomes out of range [{}, {}]".format(LOWER_LIMIT, UPPER_LIMIT),
        "description": "Autosomes whose ratio of own median to autosomal median"
        + " is out of valid range [{}, {}]".format(LOWER_LIMIT, UPPER_LIMIT),
        "scale": "Accent",
        "colour": "0, 255, 0",
    },
}


# The TABLE_CONFIG defines configs for the whole table in own section.
TABLE_CONFIG = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "id": "dragen_ploidy_estimation_metrics_table",  # ID used for the table
    "table_title": "Ploidy estimation",  # Title of the table. Used in the column config modal
    "save_file": False,  # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sortRows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}

'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The 'point.category' is used instead of 'point.x' in the 'tt_label' config, because the
Highcharts stores the keys of each sample of the given input data in the 'category' attribute
of each point object. The 'x' attribute can be seen as simple enumeration of sample's keys.
If this explanation sounds strange, here is an example of how it works in our case:
We want to plot the "1, 2, 3, ... , 22, X, Y median / Autosomal median" metrics in a linegraph.
Hovering over a single point of a line shows a description, which is defined in the 'tt_label'.
So, if we use the 'point.x' in the first curly braces, the following values will be inserted:
0, 1, 2, 3, ... , 21, 22, 23
But if we use the 'point.category' this preferable sequence is used instead:
1, 2, 3, 4, ... , 22, X, Y

More info here: https://api.highcharts.com/class-reference/Highcharts.Point#category

More configs: https://multiqc.info/docs/development/plots/#line-graphs
'''
LINEGRAPH_CONFIG = {
    "id": "dragen_ploidy_estimation_linegraph_with_ratios",
    "title": "Dragen: Chromosome median / Autosomal median ratios",
    "ylab": "Chromosome median / Autosomal median",
    "xlab": "Chromosome",
    "categories": True,  # Set to True to use x values as categories instead of numbers.
    "tt_label": "{point.category} median / Autosomal median: {point.y:.2f}x",  # Use to customise tooltip label
}


class DragenPloidyEstimationMetrics(BaseMultiqcModule):
    """Public members of the DragenPloidyEstimationMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined."""

    def add_ploidy_estimation_metrics(self):
        """The main function of the dragen ploidy metrics module.
        Public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns keys with sample names
        or an empty set if no samples were found or all samples were ignored."""

        ploidy_data = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        # Metric IDs with original consistent metric strings.
        all_metrics = {}

        for f in self.find_log_files("dragen/ploidy_estimation_metrics"):
            out = ploidy_parser(f)
            if out["success"]:
                s_name = f["s_name"]
                if s_name in ploidy_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                self.add_data_source(f, section="stats")
                all_samples[s_name].append(f)
                ploidy_data[s_name] = out["data"]  # Add/overwrite the sample.
                all_metrics.update(out["metric_IDs"])

        # Filter to strip out ignored sample names.
        ploidy_data = self.ignore_samples(ploidy_data)
        if not ploidy_data:
            return set()

        check_duplicate_samples(all_samples, log, "dragen/ploidy_estimation_metrics")

        ploidy_headers = make_headers(all_metrics)

        txt_out_data = flatten_4D_data(
            ploidy_data,
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly metric IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        self.write_data_file(txt_out_data, "dragen_ploidy")  # Write data to file.

        add_own_metrics(ploidy_data, ploidy_headers)

        ploidy_data_table = make_data_for_table(ploidy_data)  # Prepare data for tables.
        table_data, table_headers = flatten_4D_data(
            ploidy_data_table,
            headers=ploidy_headers,  # Adjust headers' metrics to data's structure.
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        gen_headers = make_gen_stats_headers(table_headers)
        # If gen_headers is empty, then all metrics from data will be included in table.
        if gen_headers:
            self.general_stats_addcols(table_data, gen_headers, namespace=NAMESPACE)

        # Prepare html table.
        own_table_headers = make_own_headers(table_headers)
        self.add_section(
            name="Ploidy Estimator",
            anchor="dragen-ploidy-estimation-metrics-table",
            description="The Ploidy Estimator results, including each normalized "
            "per-contig median coverage. Press the `Help` button for details.",
            helptext="The Ploidy Estimator runs by default. The Ploidy Estimator "
            "uses reads from the mapper/aligner to calculate the sequencing depth of "
            "coverage for each autosome and allosome in the human genome. The sex "
            "karyotype of the sample is then estimated using the ratios of the median "
            "sex chromosome coverages to the median autosomal coverage. The sex karyotype "
            "is estimated based on the range the ratios fall in. If the ratios are outside "
            "all expected ranges ,then the Ploidy Estimator does not determine a sex karyotype. "
            "Ploidy estimation can fail if the type of input sequencing data cannot be "
            "determined or if there is not sufficient sequencing coverage in the autosomes. "
            "When ploidy estimation fails the estimated median coverage values will be zero. "
            "When both tumor and matched normal reads are provided as input, the Ploidy Estimator "
            "only estimates sequencing coverage and sex karyotype for the matched normal sample "
            "and ignores the tumor reads. If only tumor reads are provided as input, the Ploidy "
            "Estimator estimates sequencing coverage and sex karyotype for the tumor sample.",
            plot=table.plot(
                table_data, own_table_headers, {config: val for config, val in TABLE_CONFIG.items() if val is not None}
            ),
        )
        # Add own section with linegraph for "Chromosome median / Autosomal median" metrics.
        linegraph_section = make_data_for_linegraph_section(ploidy_data)
        if linegraph_section["data"]:  # Append only if not empty.
            self.add_section(
                name="Chromosome median / Autosomal median",
                anchor="dragen-ploidy-estimation-metrics-linegraph",
                description="'Autosome|Allosome median / Autosomal median' ratios from .ploidy_estimation_metrics.csv files. "
                + linegraph_section["description"],
                plot=linegraph.plot(linegraph_section["data"], LINEGRAPH_CONFIG),
            )
        # Report found info/warnings/errors. You can disable it anytime, if it is not wanted.
        make_log_report(log_data, log, "ploidy_estimation_metrics")

        return ploidy_data.keys()


def make_data_for_table(DATA):
    """Modify parsed data before storing in tables."""

    # Inner samples can be concatenated with file samples
    # and replaced by empty string in the second column.

    # Sections can be analyzed and restructured.

    # Use INCLUDED_SECTIONS as inclusion checklist.
    return {
        sample: {
            section: {
                region: {metric: value for metric, value in data.items()}
                for region, data in regions.items()
                if (ANY in INCLUDED_SECTIONS and ANY in INCLUDED_SECTIONS[ANY])
                or (
                    section in INCLUDED_SECTIONS
                    and (region in INCLUDED_SECTIONS[section] or ANY in INCLUDED_SECTIONS[section])
                )
            }
            for section, regions in sections.items()
            if ANY in INCLUDED_SECTIONS or section in INCLUDED_SECTIONS
        }
        for sample, sections in DATA.items()
    }


def add_own_metrics(ploidy_data, ploidy_headers):
    """Inserts calculated data into ploidy_data.
    Appends new "metrics" to the ploidy_headers."""

    metric_id = "autosome median to autosomal median ratios out of range"
    header = {config: value for config, value in SINGLE_HEADER.items() if value is not None}
    header.update(EXTRA_METRICS[metric_id])
    ploidy_headers[metric_id] = header

    # For the "1,2,3,...,22 median / autosomal median" metrics.
    AUTOSOME_RGX = re.compile("^(\d+) median / autosomal median$")

    for sample in ploidy_data:
        if PLOIDY_ESTIMATION in ploidy_data[sample] and EMPTY in ploidy_data[sample][PLOIDY_ESTIMATION]:
            metrics = ploidy_data[sample][PLOIDY_ESTIMATION][EMPTY]
            value = ""
            for metric in metrics:
                autosome_match = AUTOSOME_RGX.search(metric)
                if autosome_match:
                    aut_val = metrics[metric]
                    # Though checking the metric's value might seem redundant, it is necessary to guarantee
                    # that MultiQC will not crash if some string (eg NA) is present instead of a number.
                    if isinstance(aut_val, (float, int)) and (aut_val < LOWER_LIMIT or aut_val > UPPER_LIMIT):
                        value += autosome_match.group(1) + "; "
            if value:
                # Delete the last space but not the semicolon. Without a semicolon at the
                # end of the string a single found autosome would be converted to float.
                # There is no need to check/create a key for ADDITIONAL_SECTION, because
                # each sample references a defaultdict(lambda: defaultdict(dict)).
                ploidy_data[sample][ADDITIONAL_SECTION][EMPTY][metric_id] = value[:-1]
            else:
                ploidy_data[sample][ADDITIONAL_SECTION][EMPTY][metric_id] = "None"


def make_data_for_linegraph_section(ploidy_data):
    """Create data with '1,2,3,...,22 median / Autosomal median' metrics for linegraph."""

    # OrderedDict to guarantee the order of metrics. Simple dict might cause wrong output
    # if python implementation does not guarantee the preservation of order (eg python 3.6).
    data = defaultdict(OrderedDict)

    # Tell users how absent metrics are handled.
    description = ""

    # The CHROMS is a list of pairs with "1,2,...,22,x,y median / autosomal median" metrics and
    # associated chromosomes. Must be in lower case, because of how make_consistent_metric works.
    # The sole purpose of this array is to guarantee the sequence of metrics. Creating the output data
    # by iterating through the input data would result in wrong output if the first catched sample has
    # a sequence of metrics, which is distinct from the sequence of the same metrics in other samples.
    CHROMS = [
        (str(chrom) + " median / autosomal median", str(chrom).upper()) for chrom in list(range(1, 23)) + ["x", "y"]
    ]
    for sample in ploidy_data:
        if PLOIDY_ESTIMATION in ploidy_data[sample] and EMPTY in ploidy_data[sample][PLOIDY_ESTIMATION]:
            metrics = ploidy_data[sample][PLOIDY_ESTIMATION][EMPTY]
            for metric, chrom in CHROMS:
                # Check if metric is in the sample and its value is a number. MultiQC can crash without the latter.
                if metric in metrics and isinstance(metrics[metric], (float, int)):
                    data[sample][chrom] = metrics[metric]
                # Additional safety measure to prevent producing broken output
                # in case the first catched sample does not contain all metrics.
                else:
                    data[sample][chrom] = -1
                    description = (
                        "Ratio of -1 means that metric is not present in the sample or its value is not a number."
                    )
    return {"data": data, "description": description}


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-file-prefix>.ploidy_estimation_metrics.csv

debug/info:
- invalid_file_lines, which do not conform to:
  <section>,<RG/Sample/Empty>,<metric>,<value>
- unknown_sections are not present in the SUPPORTED_SECTIONS.
- unknown_second_columns are not present in the SUPPORTED_SECTIONS.
- unknown_metrics are those, which are not part of the STD_METRICS.
  Their headers are incomplete, hence uninformative and ugly table's columns.
- unusual_values are those except for int/float and a combo of X/Y/0.
'''
log_data = {
    "invalid_file_names": defaultdict(list),
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_metrics": [],
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
    "unknown_second_columns": defaultdict(set),
    "unknown_sections": set(),
}


def make_headers(metric_IDs):
    """Build headers from metric_IDs."""
    headers = {}

    # The AUTOSOME_RGX is for the "1,2,3,...,22 median / autosomal median" metrics.
    AUTOSOME_RGX = re.compile("^(\d+) median / autosomal median$")

    for metric_id in metric_IDs:
        original_metric = metric_IDs[metric_id]
        configs = {config: val for config, val in SINGLE_HEADER.items() if val is not None}

        median_ratio_match = AUTOSOME_RGX.search(metric_id)
        if median_ratio_match:
            autosome_n = median_ratio_match.group(1)
            configs["order_priority"] = RESERVED + int(autosome_n)  # Order ascendingly.
            configs["title"] = autosome_n + "Med/AutMed"
            configs["description"] = (
                "Median coverage of autosome " + autosome_n + " devided by the autosomal median coverage."
            )
            _configs = STD_METRICS["autosome_n median / autosomal median"]
            configs.update(
                {config: val for config, val in _configs.items() if config in SINGLE_HEADER or config in EXTRA_HEADER}
            )
            if autosome_n in _configs["extra"] and _configs["extra"][autosome_n]:
                configs.update(_configs["extra"][autosome_n])
        else:
            if metric_id in STD_METRICS:
                configs.update(STD_METRICS[metric_id])
            else:
                log_data["unknown_metrics"].append(original_metric)
                configs["title"] = original_metric
                configs["description"] = original_metric

        headers[metric_id] = configs

    return headers


def create_parser():
    """Isolation for the parsing codeblock. Returns parse_ploidy_estimation_metrics_file closure."""

    def make_consistent_metric(metric):
        """Tries to guarantee consistency to a certain degree.
        Metric-specific peculiarities can be handled here."""
        # metric will be stored in .txt output.
        metric = re.sub("\s+", " ", metric).strip()
        # metric_id will be used for internal storage.
        metric_id = metric.lower()
        return metric, metric_id

    def make_consistent_section(section):
        """Improve consistency of section's string."""
        return re.sub("\s+", " ", section).strip().upper()

    # The official/accepted structure of the input file:
    # <output-file-prefix>.ploidy_estimation_metrics.csv
    FILE_RGX = re.compile("(.+)\.ploidy_estimation_metrics\.csv$")

    # There is no officially defined line structure.
    # The following regex is used to check input lines.
    LINE_RGX = re.compile("^([^,]+),([^,]*),([^,]+),([^,]+)$")

    def parse_ploidy_estimation_metrics_file(file_handler):
        """
        Parser for *.ploidy_estimation_metrics.csv files.

        Input:  file_handler with necessary info - file name/content/root.

        Output: {"success": False} if file name test failed. Otherwise:
                {"success": False/True,
                 "sample_name": <output-file-prefix>,
                 "data": {section: {RG/Sample: {metric: value}}},
                 "metric_IDs": {metric_ID: original_metric}
                }, with success False if all file's lines are invalid.

        File example:

        T_SRR7890936_50pc.ploidy_estimation_metrics.csv
        PLOIDY ESTIMATION,,Autosomal median coverage,55.63
        PLOIDY ESTIMATION,,X median coverage,27.44
        PLOIDY ESTIMATION,,Y median coverage,0.00
        PLOIDY ESTIMATION,,X median / Autosomal median,0.49
        PLOIDY ESTIMATION,,Y median / Autosomal median,0.00
        PLOIDY ESTIMATION,,Ploidy estimation,X0
        """
        file, root = file_handler["fn"], file_handler["root"]

        file_match = FILE_RGX.search(file)
        if not file_match:
            log_data["invalid_file_names"][root].append(file)
            return {"success": False}

        sample = file_match.group(1)

        success = False
        data = defaultdict(lambda: defaultdict(dict))
        metric_IDs = {}

        for line in file_handler["f"].splitlines():
            # Check the general line structure.
            line_match = LINE_RGX.search(line)

            # If line is fine then extract the necessary fields.
            if line_match:
                success = True
                section, region_sample, metric, value = line_match.groups()

            # Otherwise check if line is empty. If not then report it. Go to the next line.
            else:
                if not re.search("^\s*$", line):
                    log_data["invalid_file_lines"][root][file].append(line)
                continue

            # Modify extracted parts to improve consistency.
            consistent_section = make_consistent_section(section)
            consistent_metric, metric_id = make_consistent_metric(metric)

            # Check support for sections/regions.
            if not (ANY in SUPPORTED_SECTIONS or consistent_section in SUPPORTED_SECTIONS):
                log_data["unknown_sections"].add(section)
            if not (
                ANY in SUPPORTED_SECTIONS
                and ANY in SUPPORTED_SECTIONS[ANY]
                or (
                    consistent_section in SUPPORTED_SECTIONS
                    and (
                        region_sample in SUPPORTED_SECTIONS[consistent_section]
                        or ANY in SUPPORTED_SECTIONS[consistent_section]
                    )
                )
            ):
                log_data["unknown_second_columns"][section].add(region_sample)

            # Check the extracted value. Float in most cases. And a single karyotype.
            try:
                value = float(value)
            except ValueError:
                # Int is also checked, though is was not found in real data.
                try:
                    value = int(value)
                except ValueError:
                    if not re.search("^[XY0]+$", value.strip(), re.IGNORECASE):
                        log_data["unusual_values"][root][file][metric] = value

            metric_IDs[metric_id] = consistent_metric
            data[consistent_section][region_sample][metric_id] = value

        # Return results of parsing the input file.
        return {
            "success": success,
            "sample_name": sample,
            "data": data,
            "metric_IDs": metric_IDs,
        }

    # Return the parser closure.
    return parse_ploidy_estimation_metrics_file


ploidy_parser = create_parser()
