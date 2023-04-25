'''"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
Below is a set of common tools, which might help to write DRAGEN QC Metrics modules.

This comment describes some construction components of a single dragen module.
It also includes noticed common patterns, observations.


Officially defined QC metrics output table format:
Section    RG/Sample    Metric    Count/Ration/Time    Percentage/Seconds

Source: https://support-docs.illumina.com/SW/DRAGEN_v41/Content/SW/DRAGEN/QCMetricsCoverageReports.htm

Real files' formats can vary:

- Sections have variations in composition and staticness.
  * Can be constant: coverage_metrics with 'COVERAGE SUMMARY'
  * Or quite variable:
    * subsections/prefixes/suffixes (eg vc_metrics, mapping_metrics).
    * phenotypes as extra prefixes: 'TUMOR', 'NORMAL'
  * Some dont have a section:
    Eg overall_mean_cov_metrics has the 'Average alignment coverage over <source>' metric in the 1st column.
    But the 'TUMOR' and 'NORMAL' prefixes can be also appended to that metric.

- RG/Sample can be either:
  * DRAGEN-defined: empty or some static region (eg WGS).
  * User-defined: sample names.
  * Non existent (eg overall_mean_cov_metrics)

- Metrics are also quite variant:
  * Substrings:
    * Varying regions/values/metainfo (eg mapping_metrics, coverage_metrics).
    * Space padding within metrics (eg coverage_metrics).
  * DRAGEN version dependency as stated in mapping_metrics.py:
    "Dragen 3.9 has replaced 'rRNA filtered reads' with 'Adjustment of reads matching filter contigs'"
  * User-defined sample names instead of metrics (eg cnv_metrics).
  * Probably smth else.

- Values are mostly int and float. NA, inf or karyotype(eg XX/XY/X0) are also present.
  Extra padded space can be included (eg overall_mean_cov_metrics)


DRAGEN QC Metrics Module includes:

Main components:

  - Main class:
    interface to features provided by the BaseMultiqcModule and own functionality.
    Methods and properties of DRAGEN modules defined in MultiqcModule(dragen.py) are also
    available from the main dragen class through self. Unfortunately it is not possible
    to completely isolate dragen modules from each other, so some cautiousness is needed:
    * main class can overwrite public members of main classes from other modules.
    * try to create as little methods/properties as possible. Especially those which can
      consume a lot of space, because references to such members will still exist after
      single module's execution has finished.
    * use only explicitly defined internal dragen interface(eg getter methods) to quiry data.

  - CSV Parser:
    extracts data from specific CSV-file with metrics.
    The structure of collected data depends on the file's format.

  - Configurations:
    Highcharts plots: headers for tables, configs for whole tables ...

Optional components:

  - Data transformator:
    tools (eg util functions, own methods/global funcs) to restructure parsed data.
    It might happen that a single file contains differently structured rows.
    Internal data restructuring could be necessary in order to properly and nicely
    present metrics in html table.
    - Data flattener:
      Transforms data from higher dimension to 2D - {sample: {metric: value}}
      May be useful in writing a prototype.
    - Transformation for storing in txt output.
    - Transformation for tables/linegraphs/etc

  - Headers creation:
    the only plot which is available by default and can represent data completely is the table.
    In some cases metrics/files can have varying and rather complicated structure. It might be
    cumbersome to manually define configs in such cases, so additional functionality could help.

  - Logging:
    Collect and report found debug/info/warnings/errors according to some predefined schema.
    Direct reporting can be also considered as a "component", a distributed one, sort of.
    It is much easier than the former and reports with own module name, but can produce chaotic log.

  - JSON Parser:
    JSON data can be collected separately from CSV. After that it can be used to improve html output.

  - Data checker:
    It might happen that collected data needs to be thoroughly examined to guarantee data quality.
    This component may be useful if many specific data aspects have to be checked.
    Currently it is part of parser.

This structuring is just an attempt to organize and simplify the whole dragen module, which
could improve readability and make modules more extensible/maintainable in the long run.
'''

import re
from collections import OrderedDict, defaultdict
from multiqc import config


read_format = "{:,.1f}"
if config.read_count_multiplier == 1:
    read_format = "{:,.0f}"
# read_format += '&nbsp;' + config.read_count_prefix

base_format = "{:,.1f}&nbsp;"
if config.base_count_multiplier == 1:
    base_format = "{:,.0f}"
elif config.base_count_multiplier == 0.000000001:
    base_format = "{:,.2f}"
# base_format += '&nbsp;' + config.base_count_prefix


"""______________________    Headers    ______________________"""


class Metric:
    def __init__(
        self,
        id,
        title,
        in_genstats=None,
        in_own_tabl=None,
        unit=None,
        descr=None,
        fmt=None,
        modify=None,
        namespace=None,
        precision=None,
        the_higher_the_worse=False,
    ):
        self.id = id
        self.title = title
        self.in_genstats = in_genstats
        self.in_own_tabl = in_own_tabl
        self.unit = unit
        self.descr = descr
        self.fmt = fmt
        self.modify = modify
        self.namespace = namespace
        self.precision = precision
        self.the_higher_the_worse = the_higher_the_worse


def make_headers(parsed_metric_ids, metrics):
    # Init general stats table
    genstats_headers = OrderedDict()

    # Init headers for an own separate table
    own_tabl_headers = OrderedDict()

    for metric in metrics:
        col = dict(
            title=metric.title,
            description=metric.descr,
            min=0,
        )

        # choosing color based on metric namespace, guessing namespace by unit
        if metric.unit == "reads":
            col["scale"] = "Greens"
        elif metric.unit == "bases":
            col["scale"] = "Greys"
        elif metric.unit == "bp":
            col["scale"] = "Purples"
        elif metric.unit == "x" or metric.namespace == "Coverage":
            col["scale"] = "Blues"
            if metric.unit == "%":
                col["scale"] = "RdBu"
                if metric.the_higher_the_worse:
                    col["scale"] = "YlOrRd"
        elif metric.namespace == "Variants":
            col["scale"] = "Purples"

        if metric.id + " pct" in parsed_metric_ids:
            # if % value is available, showing it instead of the number value; the number value will be hidden
            pct_col = dict(
                col,
                description=metric.descr.replace(", {}", "").replace("Number of ", "% of "),
                max=100,
                suffix="%",
            )
            if metric.unit == "reads":
                pct_col["scale"] = "RdYlGn"
            elif metric.unit == "bases":
                pct_col["scale"] = "RdGy"
            elif metric.unit == "x" or metric.namespace == "Coverage":
                pct_col["scale"] = "RdBu"
            if metric.the_higher_the_worse:
                pct_col["scale"] = "YlOrRd"

            if metric.precision is not None:
                pct_col["format"] = "{:,." + str(metric.precision) + "f}"

            if metric.in_own_tabl is not None:
                show = metric.in_own_tabl == "%"
                own_tabl_headers[metric.id + " pct"] = dict(pct_col, hidden=not show)
            if metric.in_genstats is not None and metric.in_own_tabl != "#":
                genstats_headers[metric.id + " pct"] = dict(pct_col, hidden=metric.in_genstats == "hid")

        if metric.unit == "reads":
            col["description"] = col["description"].format(config.read_count_desc)
            col["title"] = config.read_count_prefix + " " + col["title"]
            col["modify"] = lambda x: x * config.read_count_multiplier
            col["shared_key"] = "read_count"
            col["format"] = read_format
        elif metric.unit == "bases":
            col["description"] = col["description"].format(config.base_count_desc)
            col["title"] = config.base_count_prefix + " " + col["title"]
            col["modify"] = lambda x: x * config.base_count_multiplier
            col["shared_key"] = "base_count"
            col["format"] = base_format
        elif metric.unit == "len":
            col["suffix"] = " bp"
            col["format"] = "{:,.0f}"
        elif metric.unit == "x":
            col["suffix"] = " x"
            col["format"] = "{:,.1f}"
        elif metric.unit == "%":
            col["suffix"] = " %"
            col["format"] = "{:,.1f}"
        elif any(
            (metric.descr and metric.descr.startswith(pref))
            for pref in ("Total number of ", "The number of ", "Number of ")
        ):
            col["format"] = "{:,.0f}"
        if metric.precision is not None:
            col["format"] = "{:,." + str(metric.precision) + "f}"

        if metric.in_own_tabl is not None:
            own_tabl_headers[metric.id] = dict(col, hidden=metric.in_own_tabl != "#")
        if metric.in_genstats is not None:
            if metric.in_genstats == "#":
                genstats_headers[metric.id] = dict(col)
            elif metric.in_genstats == "hid" and metric.in_own_tabl == "#":
                genstats_headers[metric.id] = dict(col, hidden=metric.in_genstats == "hid")

    return genstats_headers, own_tabl_headers


def order_headers(headers):
    """Produces a shallow copy with ordered single headers from the input.
    Or returns the input if no header contains the 'order_priority' config."""
    """
    Collected "order_priority" values are stored as groups to resolve
    conflicts, which may arise if several metrics have the same priority.
    The order of headers within the same priority group is not specified.
    """
    indexes = set()
    ordered_headers = defaultdict(list)
    not_ordered_headers = {}
    for metric in headers:
        configs = headers[metric].copy()
        if "order_priority" in headers[metric] and isinstance(headers[metric]["order_priority"], (int, float)):
            order_priority = headers[metric]["order_priority"]
            indexes.add(order_priority)
            ordered_headers[order_priority].append((metric, configs))
        else:
            not_ordered_headers[metric] = configs

    # Convert to list and sort in ascending order.
    if indexes:
        indexes = list(indexes)
        indexes.sort()
    # Otherwise there is nothing to order. Return the input.
    else:
        return headers

    output_headers = OrderedDict()
    for index in indexes:
        for metric, configs in ordered_headers[index]:
            output_headers[metric] = configs

    output_headers.update(not_ordered_headers)
    return output_headers


# STD_TABLE_CONFIGS contains all standard table configurations:
STD_TABLE_CONFIGS = [
    "namespace",
    "title",
    "description",
    "max",
    "min",
    "ceiling",
    "floor",
    "minRange",
    "scale",
    "bgcols",
    "colour",
    "suffix",
    "format",
    "cond_formatting_rules",
    "cond_formatting_colours",
    "shared_key",
    "modify",
    "hidden",
    "bars_zero_centrepoint",
]


def clean_headers(headers):
    """Include configurations only from STD_TABLE_CONFIGS. Extension-configurations will be excluded."""
    cleaned_headers = OrderedDict()
    for metric in headers:
        cleaned_headers[metric] = {
            config: val for config, val in headers[metric].items() if config in STD_TABLE_CONFIGS
        }
    return cleaned_headers


def make_gen_stats_headers(headers):
    gen_headers = {
        metric: configs.copy()  # Copy to be able to modify them.
        for metric, configs in headers.items()
        if metric in
        # Header is included only if "exclude" is not present or False/False-equivalent.
        [
            metric
            for metric in headers
            if not ("exclude" in headers[metric] and headers[metric]["exclude"])
        ]
    }
    return clean_headers(order_headers(gen_headers))


def make_own_headers(headers):
    """Creates table's headers for non-general section."""
    own_headers = {
        metric: configs.copy()  # Copy to be able to modify them.
        for metric, configs in headers.items()
        if metric in
        # Header is included only if "exclude_own" is not present or False/False-equivalent.
        [
            metric
            for metric in headers
            if not ("exclude_own" in headers[metric] and headers[metric]["exclude_own"])
        ]
    }
    for metric in own_headers:
        # "hidden_own" is an extension and must be stored in "hidden"
        if "hidden_own" in own_headers[metric]:
            own_headers[metric]["hidden"] = own_headers[metric]["hidden_own"]

    unhide_single_metric(own_headers)
    return clean_headers(order_headers(own_headers))


"""______________________    Data transformation    ______________________"""


def flatten_3D_data(data, **extra):
    """
    Reduces amount of dimensions in data to 2D by combining layers.

    Input:
      * data:
        dictionary with parsed metrics from coverage reports in form:
        {sample: {section: {metric_ID: value}}}

      * extra:
        includes additional parameters:
        - hide_single_section: boolean to hide single common section.
        - metric_IDs: dictionary with parsed metrics in form {metric_ID: original_metric}
        - headers: dictionary in form {metric_ID: {configuration: value}}

    Output:
      * 2D data if key parameter headers is not provided.
      * 2D data and headers if key parameter headers is provided.
    """
    out_data = defaultdict(dict)
    headers = extra["headers"] if "headers" in extra else None
    out_headers = {} if headers is not None else None

    metric_IDs = extra["metric_IDs"] if "metric_IDs" in extra else None
    hide_single_section = extra["hide_single_section"] if "hide_single_section" in extra else False
    sections = len({section for sample in data for section in data[sample]})

    for sample in data:
        for section in data[sample]:
            for metric in data[sample][section]:
                original_metric = metric
                # Try to extract the original metric.
                if metric_IDs and metric in metric_IDs:
                    original_metric = metric_IDs[metric]
                # Create new metric-ID.
                if hide_single_section and sections == 1:
                    new_id = original_metric
                else:
                    # If section is not empty.
                    if section:
                        new_id = section + " " + original_metric
                    else:
                        new_id = original_metric
                out_data[sample][new_id] = data[sample][section][metric]

                if out_headers is None:
                    continue

                if metric in headers:
                    out_headers[new_id] = headers[metric]

    if out_headers is not None:
        return out_data, out_headers
    else:
        return out_data


def flatten_4D_data(data, **extra):
    """
    Reduces amount of dimensions in data to 2D by combining layers.

    Input:
      * data:
        dictionary with parsed metrics from coverage reports in form:
        {sample: {section: {region: {metric_ID: value}}}}

      * extra:
        includes additional parameters:
        - hide_single_section: boolean to hide single common section.
        - metric_IDs: dictionary with parsed metrics in form {metric_ID: original_metric}
        - headers: dictionary in form {metric_ID: {configuration: value}}

    Output:
      * 2D data if key parameter headers is not provided.
      * 2D data and headers if key parameter headers is provided.
    """
    out_data = defaultdict(dict)
    headers = extra["headers"] if "headers" in extra else None
    out_headers = {} if headers is not None else None

    metric_IDs = extra["metric_IDs"] if "metric_IDs" in extra else None
    hide_single_section = extra["hide_single_section"] if "hide_single_section" in extra else False
    sections = len({section for sample in data for section in data[sample]})

    for sample in data:
        for section in data[sample]:
            for region in data[sample][section]:
                for metric in data[sample][section][region]:
                    original_metric = metric
                    # Try to extract the original metric.
                    if metric_IDs and metric in metric_IDs:
                        original_metric = metric_IDs[metric]
                    # Create new metric-ID.
                    if hide_single_section and sections == 1:
                        # If region is not empty.
                        if region:
                            new_id = region + " " + original_metric
                        else:
                            new_id = original_metric
                    else:
                        # Section and region are not empty.
                        if section and region:
                            new_id = section + " " + region + " " + original_metric
                        # Region is empty.
                        elif section and not region:
                            new_id = section + " " + original_metric
                        # Section is empty.
                        elif not section and region:
                            new_id = region + " " + original_metric
                        # Both are empty.
                        else:
                            new_id = original_metric
                    out_data[sample][new_id] = data[sample][section][region][metric]

                    if out_headers is None:
                        continue

                    if metric in headers:
                        out_headers[new_id] = headers[metric]

    if out_headers is not None:
        return out_data, out_headers
    else:
        return out_data


"""______________________    Logging component    ______________________"""


def check_duplicate_samples(sample_names, logger, module):
    """Check samples for duplicate names. Warn about found ones."""
    message = ""
    line1 = "\n  {} was built from the following samples:"
    line2 = "\n    {} in {}"
    for clean_sample_name in sample_names:
        if len(sample_names[clean_sample_name]) > 1:
            message += line1.format(clean_sample_name)
            for file_handler in sample_names[clean_sample_name]:
                message += line2.format(file_handler["fn"], file_handler["root"])
    if message:
        # Draw user's attention with error.
        logger.error(module + ": duplicates found.")
        # But print full info only in report.
        message = "\n\nDuplicate sample names were found. The last one overwrites previous data." + message
        logger.debug(message)


# Used to define module-specific texts for the make_parsing_log_report.
# Comments are used to separate keys. Though empty lines look better, they are removed by the "Black" formatter.
DRAGEN_MODULE_TEXTS = {
    "invalid_file_names": {
        "coverage_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output-prefix>.<coverage-region-prefix>_coverage_metrics<arbitrary-suffix>.csv\n\n"
        "The following files are not valid:\n",
        # ------------- #
        "overall_mean_cov_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output-prefix>_overall_mean_cov<arbitrary-suffix>.csv\n\n"
        "The following files are not valid:\n",
        # ------------- #
        "cnv_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output prefix>.cnv_metrics.csv\n\n"
        "The following files are not valid:\n",
        # ------------- #
        "ploidy_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output-file-prefix>.ploidy_estimation_metrics.csv\n"
        "The following files are not valid:\n",
        # ------------- #
        "roh_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output-file-prefix>.roh_metrics.csv\n"
        "The following files are not valid:\n",
        # ------------- #
        "sv_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output-file-prefix>.sv_metrics.csv\n"
        "The following files are not valid:\n",
        # ------------- #
        "vc_hethom_ratio_metrics": "\n\nThe file names must conform to the following structure:\n"
        "<output-file-prefix>.vc_hethom_ratio_metrics.csv\n"
        "The following files are not valid:\n",
    },
    "invalid_file_lines": {
        "coverage_metrics": "\n\nThe lines in files must be:\n"
        "COVERAGE SUMMARY,,<metric>,<value1> or "
        "COVERAGE SUMMARY,,<metric>,<value1>,<value2>\n"
        "The following files contain invalid lines:\n",
        # ------------- #
        "cnv_metrics": "\n\nThe lines in files must be:\n"
        "<section>,<RG/Sample/Empty>,<sample>,<value> or "
        "<section>,<RG/Sample/Empty>,<sample>,<value1>,<value2>\n"
        "The following files contain invalid lines:\n",
        # ------------- #
        "ploidy_metrics": "\n\nThe lines in files must be:\n"
        "<section>,<RG/Sample/Empty>,<metric>,<value>\n"
        "The following files contain invalid lines:\n",
        # ------------- #
        "roh_metrics": "\n\nThe lines in files must be:\n"
        "<section>,<RG/Sample/Empty>,<metric>,<value>\n"
        "The following files contain invalid lines:\n",
        # ------------- #
        "sv_metrics": "\n\nThe lines in files must be:\n"
        "<section>,<RG/Sample/Empty>,<metric>,<value> or "
        "<section>,<RG/Sample/Empty>,<metric>,<value1>,<value2>\n"
        "The following files contain invalid lines:\n",
        # ------------- #
        "vc_hethom_ratio_metrics": "\n\nThe lines in files must be:\n"
        "<section>,<contig>,<metric>,<value>\n"
        "The following files contain invalid lines:\n",
    },
    "unknown_sections": {},
    "unknown_second_columns": {},
    "unknown_metrics": {
        "coverage_metrics": "\n\nThe following metrics could not be recognized. "
        "Their headers might be incomplete, hence uninformative and ugly table's columns.\n",
        # ------------- #
        "overall_mean_cov_metrics": "\n\nOnly the 'Average alignment coverage over <file>' metric is supported.\n"
        "The following metrics could not be recognized and are not included in the report.\n",
        # ------------- #
        "ploidy_metrics": "\n\nThe following metrics could not be recognized. "
        "Their headers are incomplete, hence uninformative and ugly table's columns.\n",
        # ------------- #
        "roh_metrics": "\n\nThe following metrics could not be recognized. "
        "Their headers are incomplete, hence uninformative and ugly table's columns.\n",
        # ------------- #
        "sv_metrics": "\n\nThe following metrics are not included in the report:\n",
        # ------------- #
        "vc_hethom_ratio_metrics": "\n\nThe following metrics are not included in the report:\n",
    },
    "unusual_values": {
        "coverage_metrics": "\n\nAll metrics' values except for int, float and NA are non-standard.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "overall_mean_cov_metrics": "\n\nAll metrics' values except for float are non-standard.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "cnv_metrics": "\n\nThe CNV SUMMARY section shall contain only int and float.\n"
        "The second value of the SEX GENOTYPER section shall be float.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "ploidy_metrics": "\n\nEach value shall be either int, float or a combination of X/Y/0 for the karyotype.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "roh_metrics": "\n\nEach value shall be either int or float.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "sv_metrics": "\n\nEach value shall be either int or float.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "vc_hethom_ratio_metrics": "\n\nEach value shall be either int, float or inf.\n"
        "The following files contain non-standard values:\n",
    },
}


def make_log_report(log_data, logger, module):
    """The only purpose of this function is to create a readable and
    informative log output about found debug/info/warnings/errors"""

    if "invalid_file_names" in log_data and log_data["invalid_file_names"]:
        if module in DRAGEN_MODULE_TEXTS["invalid_file_names"]:
            log_message = DRAGEN_MODULE_TEXTS["invalid_file_names"][module]
        else:
            log_message = "\n\nThe following files are not valid:\n"

        for root in log_data["invalid_file_names"]:
            log_message += "  " + root + ":\n"
            for file in log_data["invalid_file_names"][root]:
                log_message += "    " + file + "\n"

        # Draw user's attention with error.
        logger.error("dragen/" + module + ": invalid file names found.")
        # But print full info only in report.
        logger.debug(log_message + "\n")

    if "unknown_sections" in log_data and log_data["unknown_sections"]:
        if module in DRAGEN_MODULE_TEXTS["unknown_sections"]:
            log_message = DRAGEN_MODULE_TEXTS["unknown_sections"][module]
        else:
            log_message = "\n\nThe following sections could not be recognized.\n"

        for section in log_data["unknown_sections"]:
            log_message += "  " + section + "\n"

        # Draw user's attention with warning.
        logger.warning("dragen/" + module + ": unknown sections found.")
        # But print full info only in report.
        logger.debug(log_message + "\n")

    if "unknown_second_columns" in log_data and log_data["unknown_second_columns"]:
        if module in DRAGEN_MODULE_TEXTS["unknown_second_columns"]:
            log_message = DRAGEN_MODULE_TEXTS["unknown_second_columns"][module]
        else:
            log_message = "\n\nThe following RG/Samples could not be recognized.\n"

        for section in log_data["unknown_second_columns"]:
            log_message += "  " + section + ":\n"
            for rg_sample in log_data["unknown_second_columns"][section]:
                rg_sample = rg_sample if rg_sample else "Empty region"
                log_message += "    " + rg_sample + "\n"

        # Draw user's attention with warning.
        logger.warning("dragen/" + module + ": unknown RG/Samples found.")
        # But print full info only in report.
        logger.debug(log_message + "\n")

    if "invalid_file_lines" in log_data and log_data["invalid_file_lines"]:
        if module in DRAGEN_MODULE_TEXTS["invalid_file_lines"]:
            log_message = DRAGEN_MODULE_TEXTS["invalid_file_lines"][module]
        else:
            log_message = "\n\nThe following files contain invalid lines:\n"

        for root in log_data["invalid_file_lines"]:
            log_message += "  " + root + ":\n"
            for file in log_data["invalid_file_lines"][root]:
                log_message += "    " + file + ":\n"
                for line in log_data["invalid_file_lines"][root][file]:
                    log_message += "      " + line + "\n"

        # Draw user's attention with warning.
        logger.warning("dragen/" + module + ": invalid file lines found.")
        # But print full info only in report.
        logger.debug(log_message + "\n")

    if "unknown_metrics" in log_data and log_data["unknown_metrics"]:
        if module in DRAGEN_MODULE_TEXTS["unknown_metrics"]:
            log_message = DRAGEN_MODULE_TEXTS["unknown_metrics"][module]
        else:
            log_message = "\n\nThe following metrics could not be recognized:\n"

        for metric in log_data["unknown_metrics"]:
            log_message += "  " + metric + "\n"

        # Draw user's attention with warning.
        logger.warning("dragen/" + module + ": unknown metrics found.")
        # But print full info only in report.
        logger.debug(log_message + "\n")

    if "unusual_values" in log_data and log_data["unusual_values"]:
        if module in DRAGEN_MODULE_TEXTS["unusual_values"]:
            log_message = DRAGEN_MODULE_TEXTS["unusual_values"][module]
        else:
            log_message = "\n\nThe following files contain non-standard values:\n"

        for root in log_data["unusual_values"]:
            log_message += "  " + root + ":\n"
            for file in log_data["unusual_values"][root]:
                log_message += "    " + file + ":\n"
                for metric in log_data["unusual_values"][root][file]:
                    log_message += "      " + metric + " = " + log_data["unusual_values"][root][file][metric] + "\n"

        # Draw user's attention with warning.
        logger.warning("dragen/" + module + ": unusual values found.")
        # But print full info only in report.
        logger.debug(log_message + "\n")


"""______________________    Miscellaneous    ______________________"""


def exist_and_number(data, *metrics):
    return all(isinstance(data.get(m, None), int) or isinstance(data.get(m, None), float) for m in metrics)


def unhide_single_metric(headers):
    """Solve a table bug: if table has only 1 metric and this metric is hidden, then there
    is no user-friendly way to show it, because the "Configure columns" button is absent."""
    if len(headers) == 1:
        metric = list(headers.keys())[0]
        headers[metric]["hidden"] = False
