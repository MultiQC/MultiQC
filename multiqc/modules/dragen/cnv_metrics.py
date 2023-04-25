'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
This module gathers CNV metrics data and prepares it for the output report.
It relies on the following official sources:
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/CNVMetrics.htm
https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/CNVSampleCorrelationSexGenotyper.htm
'''

import logging
import re
from collections import defaultdict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file

from .utils import check_duplicate_samples, flatten_4D_data, make_gen_stats_headers, make_log_report, make_own_headers

# Initialise the logger.
log = logging.getLogger(__name__)


NAMESPACE = "DRAGEN CNV"


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The SINGLE_HEADER is used to define the common settings for all headers:
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
    "format": "{:,.1f}",  # Value format string - default 1 decimal place
    "cond_formatting_rules": None,  # Rules for conditional formatting table cell values.
    "cond_formatting_colours": None,  # Styles for conditional formatting of table cell values
    "shared_key": None,  # See the link for description
    "modify": None,  # Lambda function to modify values, special case, see below
    "hidden": True,  # Set to True to hide the column on page load
    "hidden_own": True,  # For non-general plots in the CNV section.
    "exclude": True,  # True to exclude all metrics from the general html table.
}
'''""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The EXTRA_HEADER contains the logical keys of the SINGLE_HEADER. These are just extensions,
but can be used as if they were a part of the SINGLE_HEADER. Extra configurations are not
set by default. Only the SINGLE_HEADER is used as a basis for a header. But you can always
append an extra config(eg hidden_own, exclude) to the SINGLE_HEADER to overcome this issue.

Hints and rules:
- "hidden" is used for own CNV section if "hidden_own" is not provided.

- "exclude" and "exclude_own" may be useful in cases where html size matters.
  Eg. many samples but only few metrics of interest.

- "order_priority" can be used to arrange table's columns. Only int and float are valid.
'''
EXTRA_HEADER = {
    "hidden_own": None,  # For non-general plot in the CNV section.
    "exclude": None,  # True to exclude metric from the general html table.
    "exclude_own": None,  # True to exclude from own CNV section html table.
    "order_priority": None,  # Used to specify columns' order in all tables.
}

'''"""""""""""""""""""""""""""""""""""""""""""""""
Officially defined file format could not be found.
According to examined real data, there are 5 columns in total:
Section        Void-column    Metric/Sample    Value   Pct
SEX GENOTYPER                    samples
CNV SUMMARY                      metrics

The second column does not contain anything. It is unknown whether
it's the final constant structure or it was reserved for future use.
For some reason DRAGEN violates its own official standard file format,
because it stores sample names instead of metrics in the third column.
So the module is more complicated and uglier than it could have been.

The following structure was chosen to store parsed CNV data:
{sample: {section: {region_or_sample: {metric_or_sample: value}}}}

If you wish to add new sections, take a look at the make_consistent_section.
'''

## Sections:
SEX_GENOTYPER = "SEX GENOTYPER"
CNV_SUMMARY = "CNV SUMMARY"

## RG/Sample aka second column:
EMPTY = ""
WGS = "WGS"  # Not detected/used. Just an example.

## Special case. Token for any section/region/sample.
ANY = object()

"""
Only INCLUDED_SECTIONS will be shown in html table.
Can be used as a quick but unreliable way to control html table.
"""
INCLUDED_SECTIONS = {
    SEX_GENOTYPER: [EMPTY],
    CNV_SUMMARY: [EMPTY],
}

# INCLUDED_SECTIONS = {ANY:[ANY]}  # Include all sections/regions/samples.

# Used to detect and report what is not properly supported by the module.
# Works as INCLUDED_SECTIONS.
SUPPORTED_SECTIONS = {
    SEX_GENOTYPER: [EMPTY],
    CNV_SUMMARY: [EMPTY],
}

'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""
The STD_METRICS is the container for the standard metrics.
They are designed to be seen as simple identifiers of the well known metrics.
You can set any real/virtual configuration defined in the SINGLE_HEADER.

Some rules:

- "modify" must conform to the following general structure:
  lambda x: (arbitrary expression) if isinstance(x, str) else (math expression)
  The x is either number or string. The latter is some value, which could not be
  converted to int or float. NA for instance. Such values can be modified in the
  (arbitrary expression) if needed. If not then just return the input.

- If you wish to add new metrics take a look at the make_consistent_metric.
'''

# Used to create a unique ID if metric has two values.
V2 = " pct"

RESERVED_START = 100  # Used to order panel of samples. Do not exceed it.

# Used to make a number "narrower" to not overlap with adjacent right columns.
MODIFY = lambda x: x if isinstance(x, str) else x * 0.000001
PREFIX = "M"
DESCR_UNIT = "millions"

# The following 2 metrics from the SEX_GENOTYPER section are stored as a single sample
# name in real data, so they have to get some common metric names instead.
KARYOTYPE_CASE_SAMPLE = "karyotype case sample"
CONFIDENCE_CASE_SAMPLE = "confidence case sample"

STD_METRICS = {
    KARYOTYPE_CASE_SAMPLE: {
        "order_priority": 0,
        "exclude": False,
        "hidden_own": False,
        "title": "Kar",
        "colour": "0, 255, 0",
        "hidden": False,
        "description": "Estimated sex karyotype for the case sample.",
        "cond_formatting_rules": {
            "red": [{"s_contains": ""}],  # Each value is red by default.
            "green": [{"s_eq": "XX"}, {"s_eq": "YX"}, {"s_eq": "XY"}],
        },
        "cond_formatting_colours": [
            {"red": "#FF0000"},
            {"green": "#00FF00"},
        ],
    },
    CONFIDENCE_CASE_SAMPLE: {
        "order_priority": 0.1,
        "exclude": False,
        "hidden_own": False,
        "title": "KarConf",
        "colour": "0, 255, 0",
        "hidden": False,
        "suffix": " %",
        "format": "{:,.4f}",
        "max": 1.0,
        "cond_formatting_rules": {
            "red": [{"lt": 0}, {"gt": 1}],
        },
        "cond_formatting_colours": [
            {"red": "#FF0000"},
        ],
        "description": "Confidence metric ranging from 0.0 to 1.0. If the case sample sex is specified, this metric is 0.0.",
    },
    "bases in reference genome": {
        "order_priority": 1,
        "title": config.base_count_prefix + " bases",
        "scale": "RdYlGn",
        "colour": "0, 0, 255",
        "format": "{:,.0f}",
        "description": "Bases ({}) in reference genome in use.".format(config.base_count_desc),
        "modify": lambda x: x if isinstance(x, str) else x * config.base_count_multiplier,
    },
    "average alignment coverage over genome": {
        "order_priority": 2,
        "hidden_own": False,
        "title": "AvgCov",
        "scale": "YlOrRd",
        "colour": "150, 50, 0",
        "suffix": " %",
        "description": "Average alignment coverage over genome.",
    },
    "number of alignment records": {
        "order_priority": 3,
        "hidden_own": False,
        "title": PREFIX + " Aln records",
        "scale": "BrBG",
        "colour": "0, 55, 100",
        "format": "{:,.0f}",
        "modify": MODIFY,
        "description": "Number ({}) of alignment records processed.".format(DESCR_UNIT),
    },
    "coverage uniformity": {
        "order_priority": 4,
        "hidden_own": False,
        "title": "CovUnif",
        "scale": "Reds",
        "colour": "0, 0, 255",
        "suffix": " %",
        "format": "{:,.2f}",
        "description": "Coverage uniformity.",
    },
    "number of target intervals": {
        "order_priority": 5,
        "title": "Trg intervals",
        "scale": "Greens",
        "colour": "0, 255, 200",
        "format": "{:,.0f}",
        "description": "Number of target intervals.",
    },
    "number of segments": {
        "order_priority": 6,
        "title": "Segments",
        "scale": "Blues",
        "colour": "250, 0, 0",
        "format": "{:,.0f}",
        "description": "Number of segments.",
    },
    "number of filtered records (total)": {
        "order_priority": 7.0,
        "title": PREFIX + " Total",
        "scale": "Purples",
        "colour": "255, 0, 255",
        "format": "{:,.2f}",
        "modify": MODIFY,
        "description": "Number ({}) of filtered records (total).".format(DESCR_UNIT),
    },
    "number of filtered records (total)"
    + V2: {
        "order_priority": 7.1,
        "title": "Total%",
        "scale": "Oranges",
        "colour": "255, 0, 255",
        "suffix": " %",
        "description": "Percentage of filtered records (total).",
    },
    "number of filtered records (duplicates)": {
        "order_priority": 7.2,
        "title": "Duplicates",
        "colour": "255, 0, 255",
        "scale": "Greens",
        "format": "{:,.2f}",
        "description": "Number of filtered records (due to duplicates).",
    },
    "number of filtered records (duplicates)"
    + V2: {
        "order_priority": 7.3,
        "title": "Duplicates%",
        "colour": "255, 0, 255",
        "suffix": " %",
        "description": "Percentage of filtered records (due to duplicates).",
    },
    "number of filtered records (mapq)": {
        "order_priority": 7.4,
        "title": PREFIX + " MAPQ",
        "colour": "255, 0, 255",
        "scale": "Reds",
        "format": "{:,.2f}",
        "modify": MODIFY,
        "description": "Number ({}) of filtered records (due to MAPQ).".format(DESCR_UNIT),
    },
    "number of filtered records (mapq)"
    + V2: {
        "order_priority": 7.5,
        "title": "MAPQ%",
        "colour": "255, 0, 255",
        "suffix": " %",
        "description": "Percentage of filtered records (due to MAPQ).",
    },
    "number of filtered records (unmapped)": {
        "order_priority": 7.6,
        "title": PREFIX + " Unmapped",
        "scale": "PiYG",
        "colour": "255, 0, 255",
        "format": "{:,.2f}",
        "modify": MODIFY,
        "description": "Number ({}) of filtered records (due to being unmapped).".format(DESCR_UNIT),
    },
    "number of filtered records (unmapped)"
    + V2: {
        "order_priority": 7.7,
        "title": "Unmapped%",
        "colour": "255, 0, 255",
        "suffix": " %",
        "description": "Percentage of filtered records (due to being unmapped).",
    },
    "number of amplifications": {
        "order_priority": 8.0,
        "title": "Amplifs",
        "scale": "RdYlGn",
        "colour": "255, 255, 0",
        "format": "{:,.0f}",
        "description": "Number of amplifications.",
    },
    "number of passing amplifications": {
        "order_priority": 8.1,
        "title": "Amplifs PASS",
        "scale": "PiYG",
        "colour": "255, 255, 0",
        "format": "{:,.0f}",
        "description": "Number of PASS amplifications.",
    },
    "number of passing amplifications"
    + V2: {
        "order_priority": 8.2,
        "title": "Amplifs PASS%",
        "colour": "255, 255, 0",
        "suffix": " %",
        "description": "Percentage of PASS amplifications.",
    },
    "number of deletions": {
        "order_priority": 9.0,
        "title": "Deletions",
        "scale": "Reds",
        "colour": "0, 0, 255",
        "format": "{:,.0f}",
        "description": "Number of deletions.",
    },
    "number of passing deletions": {
        "order_priority": 9.1,
        "title": "Deletions PASS",
        "scale": "BrBG",
        "colour": "0, 0, 255",
        "format": "{:,.0f}",
        "description": "Number of PASS deletions.",
    },
    "number of passing deletions"
    + V2: {
        "order_priority": 9.2,
        "title": "Deletions PASS%",
        "scale": "Oranges",
        "colour": "0, 0, 255",
        "suffix": " %",
        "description": "Percentage of PASS deletions.",
    },
    "number of de novo calls": {
        "order_priority": 10.0,
        "title": "de Novo",
        "scale": "Greens",
        "colour": "255, 0, 0",
        "format": "{:,.0f}",
        "description": "Number of de novo calls.",
    },
    "number of passing de novo calls": {
        "order_priority": 10.1,
        "title": "de Novo PASS",
        "colour": "255, 0, 0",
        "format": "{:,.0f}",
        "description": "Number of passing de novo calls.",
    },
    "number of passing de novo calls"
    + V2: {
        "order_priority": 10.2,
        "title": "de Novo PASS%",
        "scale": "Purples",
        "colour": "255, 0, 0",
        "suffix": " %",
        "description": "Percentage of passing de novo calls.",
    },
}


# The TABLE_CONFIG defines configs for the whole own table.
TABLE_CONFIG = {
    "namespace": NAMESPACE,  # Name for grouping. Prepends desc and is in Config Columns modal
    "id": "dragen_cnv_metrics_table",  # ID used for the table
    "table_title": "CNV Metrics",  # Title of the table. Used in the column config modal
    "save_file": False,  # Whether to save the table data to a file
    "raw_data_fn": None,  # File basename to use for raw data file
    "sortRows": True,  # Whether to sort rows alphabetically
    "only_defined_headers": True,  # Only show columns that are defined in the headers config
    "col1_header": None,  # The header used for the first column with sample names.
    "no_beeswarm": False,  # Force a table to always be plotted (beeswarm by default if many rows)
}


class DragenCNVMetrics(BaseMultiqcModule):
    """Public members of the DragenCNVMetrics module are available to other dragen modules.
    To avoid cluttering up the global dragen space, only one method is defined."""

    def add_cnv_metrics(self):
        """The main function of the dragen CNV metrics module.
        Public members of the BaseMultiqcModule and dragen modules defined in
        MultiqcModule are available within it. Returns keys with sample names
        or an empty set if no samples were found or all samples were ignored."""

        cnv_data = {}

        # Stores cleaned sample names with references to file handlers.
        # Used later to check for duplicates and inform about found ones.
        all_samples = defaultdict(list)

        # Metric IDs with original consistent metric strings. karyotype and confidence
        # are prestored, because these metrics are just sample names in actual data.
        all_metrics = {
            KARYOTYPE_CASE_SAMPLE: "Karyotype case sample",
            CONFIDENCE_CASE_SAMPLE: "Confidence case sample",
        }

        for file in self.find_log_files("dragen/cnv_metrics"):
            out = cnv_parser(file)
            if out["success"]:
                self.add_data_source(file, section="stats")
                cleaned_sample = self.clean_s_name(out["sample_name"], file)
                all_samples[cleaned_sample].append(file)
                cnv_data[cleaned_sample] = out["data"]  # Add/overwrite the sample.
                all_metrics.update(out["metric_IDs"])

        # Filter to strip out ignored sample names.
        cnv_data = self.ignore_samples(cnv_data)
        if not cnv_data:
            return set()

        check_duplicate_samples(all_samples, log, "dragen/cnv_metrics")

        cnv_data_modified, inner_samples_metrics = restructure_data(cnv_data)

        cnv_headers = make_headers(all_metrics, inner_samples_metrics)

        txt_out_data = flatten_4D_data(
            cnv_data_modified,
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly metric IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        self.write_data_file(txt_out_data, "dragen_cnv")  # Write data to file.

        cnv_data_table = make_data_for_table(cnv_data_modified)  # Prepare data for tables.
        table_data, table_headers = flatten_4D_data(
            cnv_data_table,
            headers=cnv_headers,  # Adjust headers' metrics to data's structure.
            metric_IDs=all_metrics,  # Original strings can be used instead of somewhat ugly IDs.
            hide_single_section=True,  # True to exclude common single section from metric ID.
        )
        gen_headers = make_gen_stats_headers(table_headers)
        # If gen_headers is empty, then all metrics will be included in table.
        if gen_headers:
            self.general_stats_addcols(table_data, gen_headers, namespace=NAMESPACE)

        # Add single own table section.
        own_table_headers = make_own_headers(table_headers)
        self.add_section(
            name="CNV Metrics",
            anchor="dragen-cnv-metrics-table",
            description="DRAGEN CNV outputs metrics in CSV format. "
            "The output follows the general convention for QC metrics reporting in DRAGEN. "
            "Press the `Help` button for details.",
            helptext="The CNV metrics are output to a file with a *.cnv_metrics.csv file extension. "
            "The following list summarizes the metrics that are output from a CNV run.\n\n"
            "Sex Genotyper Metrics:\n\n"
            "* Estimated sex karyotype for the sample, along with a confidence metric ranging from 0.0 to 1.0."
            "If the sample sex is specified, this metric is 0.0.\n\n"
            "* If using a panel of normals, all panel samples are also reported.\n\n"
            "CNV Summary Metrics:\n\n"
            "* Bases in reference genome in use.\n\n"
            "* Average alignment coverage over genome.\n\n"
            "* Number of alignment records processed.\n\n"
            "* Number of filtered records (total).\n\n"
            "* Number of filtered records (due to duplicates).\n\n"
            "* Number of filtered records (due to MAPQ).\n\n"
            "* Number of filtered records (due to being unmapped).\n\n"
            "* Number of target intervals.\n\n"
            "* Number of normal samples.\n\n"
            "* Number of segments.\n\n"
            "* Number of amplifications.\n\n"
            "* Number of deletions.\n\n"
            "* Number of PASS amplifications.\n\n"
            "* Number of PASS deletions.",
            plot=table.plot(
                table_data, own_table_headers, {config: val for config, val in TABLE_CONFIG.items() if val is not None}
            ),
        )
        # Report found info/warnings/errors. You can disable it anytime, if it is not wanted.
        make_log_report(log_data, log, "cnv_metrics")

        return cnv_data.keys()


def restructure_data(data):
    """Modify data to make it better looking in tables. The output must be 4D.
    Also produces a set with inner sample names, which will be used to create headers."""

    data_out = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    # inner_samples is the container for "metrics", which are actually inner samples.
    inner_samples = set()

    for sample in data:
        for section in data[sample]:
            # SEX_GENOTYPER has sample names instead of metrics.
            if section == SEX_GENOTYPER:
                for region in data[sample][section]:
                    # If there is only one "metric", then it is the case sample.
                    if len(data[sample][section][region]) == 1:
                        case_sample = list(data[sample][section][region].keys())[0]
                        value1, value2 = data[sample][section][region][case_sample]
                        data_out[sample][section][region][KARYOTYPE_CASE_SAMPLE] = value1
                        data_out[sample][section][region][CONFIDENCE_CASE_SAMPLE] = value2
                    # Otherwise shall be a case sample and a panel of samples.
                    else:
                        case_sample_not_catched = True  # Protection against overwriting.
                        for inner_sample in data[sample][section][region]:
                            value1, value2 = data[sample][section][region][inner_sample]  # Extract values.
                            # Check if case sample. It must be a subset of the file's sample name.
                            # Otherwise there is no way to distinguish case from panels.
                            if case_sample_not_catched and re.search(inner_sample, sample, re.IGNORECASE):
                                case_sample_not_catched = False
                                data_out[sample][section][region][KARYOTYPE_CASE_SAMPLE] = value1
                                data_out[sample][section][region][CONFIDENCE_CASE_SAMPLE] = value2
                            else:
                                inner_samples.add(inner_sample)
                                data_out[sample][section][region]["karyotype panel sample " + inner_sample] = value1
                                data_out[sample][section][region]["confidence panel sample " + inner_sample] = value2
            # Other sections are just copied.
            else:
                data_out[sample][section] = {
                    rg_sample: {metric: value for metric, value in metric_values.items()}
                    for rg_sample, metric_values in data[sample][section].items()
                }
    return data_out, inner_samples


def make_data_for_table(DATA):
    """Use INCLUDED_SECTIONS as inclusion checklist."""
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


'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
The log_data is the container for found info, warnings and errors.
Collected data is stored in groups to make the output log file more readable.

warnings:
- invalid_file_names, which do not conform to the:
    <output-file-prefix>.cnv_metrics.csv

debug/info:
- invalid_file_lines, which do not conform to:
  <section>,<RG/Sample/Empty>,<metric/sample>,<value> or
  <section>,<RG/Sample/Empty>,<metric/sample>,<value1>,<value2>
- unknown_sections are not present in the SUPPORTED_SECTIONS.
- unknown_second_columns are not present in the SUPPORTED_SECTIONS.
- unknown_metrics are not present in the STD_METRICS.
- unusual_values are those except for int, float and karyotype.
'''
log_data = {
    "invalid_file_names": defaultdict(list),
    "invalid_file_lines": defaultdict(lambda: defaultdict(list)),
    "unknown_sections": set(),
    "unknown_second_columns": defaultdict(set),
    "unknown_metrics": [],
    "unusual_values": defaultdict(lambda: defaultdict(dict)),
}


def make_headers(metric_IDs, sample_IDs):
    """Build headers from collected metrics/inner_samples."""

    # Used to order found panel samples.
    order_counter = RESERVED_START

    def make_panel_headers(sample):
        """Creates headers for both values for the given panel sample."""

        nonlocal order_counter

        configs = {
            config: val for config, val in SINGLE_HEADER.items() if val is not None
        }

        configs["order_priority"] = order_counter
        configs["title"] = "Kar " + sample
        configs["colour"] = "0, 255, 255"
        configs["description"] = (
            "Estimated sex karyotype for the panel sample " + sample
        )
        configs["cond_formatting_rules"] = {
            "red": [{"s_contains": ""}],  # Each value is red by default.
            "green": [{"s_eq": "XX"}, {"s_eq": "YX"}, {"s_eq": "XY"}],
        }
        configs["cond_formatting_colours"] = [
            {"red": "#FF0000"},
            {"green": "#00FF00"},
        ]

        order_counter += 1

        configs2 = {
            config: val for config, val in SINGLE_HEADER.items() if val is not None
        }
        configs2["order_priority"] = order_counter
        configs2["title"] = "KarConf " + sample
        configs2["colour"] = "0, 255, 255"
        configs2["description"] = (
            "Confidence metric for the panel sample "
            + sample
            + " ranging from 0.0 to 1.0. If the sample sex is specified, this metric is 0.0."
        )
        configs2["suffix"] = " %"
        configs2["format"] = "{:,.4f}"
        configs2["max"] = 1.0
        configs2["cond_formatting_rules"] = {
            "red": [{"lt": 0}, {"gt": 1}],
        }
        configs2["cond_formatting_colours"] = [
            {"red": "#FF0000"},
        ]

        order_counter += 1

        return configs, configs2

    def make_metric_header(metric, original_metric):
        """Creates a single header for the given metric."""
        configs = {
            config: val for config, val in SINGLE_HEADER.items() if val is not None
        }
        if metric in STD_METRICS:
            configs.update(STD_METRICS[metric])
        # If not, set the bare minimum.
        else:
            log_data["unknown_metrics"].append(original_metric)
            configs["title"] = original_metric
            configs["description"] = original_metric

        return configs

    headers = {}
    for metric in metric_IDs:
        headers[metric] = make_metric_header(metric, metric_IDs[metric])

    for sample in sample_IDs:
        configs, configs2 = make_panel_headers(sample)
        headers["karyotype panel sample " + sample] = configs
        headers["confidence panel sample " + sample] = configs2

    return headers


def create_parser():
    """Isolation for the parsing codeblock. Returns parse_cnv_metrics_file closure."""

    def make_consistent_metric(metric):
        """Tries to guarantee consistency to a certain degree.
        Metric-specific peculiarities can be handled here."""
        # metric will be stored in .txt output.
        metric = re.sub("\s+", " ", metric).strip()
        # metric_id will be used for internal storage.
        metric_id = metric.lower()
        return metric, metric_id

    def make_consistent_section(section):
        """Tries to guarantee consistency for section string."""
        return re.sub("\s+", " ", section).strip().upper()

    # The official/accepted structure of the input file:
    # <output-file-prefix>.cnv_metrics.csv
    FILE_RGX = re.compile("(.+)\.cnv_metrics\.csv$")

    # The following regex is used to check input lines.
    LINE_RGX = re.compile("^([^,]+),([^,]*),([^,]+),([^,]+)(,[^,]+)?$")

    def parse_cnv_metrics_file(file_handler):
        """
        Parser for *.cnv_metrics.csv files.

        Input:  file_handler with necessary info - file name/content/root.

        Output: {"success": False} if file name test failed. Otherwise:
                {"success": False/True,
                 "sample_name": <output-file-prefix>,
                 "data": {section: {RG/Sample/Empty: {metric/sample: value}}},
                 "metric_IDs": {metric_ID: original_metric}
                }, with success False if all file's lines are invalid.
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
                section, region_sample, metric, value1, value2 = line_match.groups()
                if value2:
                    value2 = value2[1:]  # Get rid of comma.

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

            # Check the extracted value. Must be int or karyotype. Float is also checked.
            try:
                value1 = int(value1)
            except ValueError:
                try:
                    value1 = float(value1)
                except ValueError:
                    if not re.search("^[XY0]+$", value1.strip(), re.IGNORECASE):
                        log_data["unusual_values"][root][file][metric] = value1

            if value2:  # Must be float. Int is also checked.
                try:
                    value2 = float(value2)
                except ValueError:
                    try:
                        value2 = int(value2)
                    except ValueError:
                        log_data["unusual_values"][root][file][metric + " (second value)"] = value2

            # SEX_GENOTYPER stores samples instead of metrics.
            if consistent_section == SEX_GENOTYPER:
                # Use original metric instead of ID to preserve letter case.
                data[consistent_section][region_sample][metric] = (value1, value2)
            else:
                metric_IDs[metric_id] = consistent_metric
                data[consistent_section][region_sample][metric_id] = value1
                if value2 is not None:
                    metric_IDs[metric_id + V2] = consistent_metric + V2
                    data[consistent_section][region_sample][metric_id + V2] = value2

        # Return results of parsing the input file.
        return {
            "success": success,
            "sample_name": sample,
            "data": data,
            "metric_IDs": metric_IDs,
        }

    # Return the parser closure.
    return parse_cnv_metrics_file


cnv_parser = create_parser()
