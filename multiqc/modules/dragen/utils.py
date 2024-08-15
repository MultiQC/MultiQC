from multiqc import config

read_format = "{:,.1f}"
if config.read_count_multiplier == 1:
    read_format = "{:,.0f}"

base_format = "{:,.1f}"
if config.base_count_multiplier == 1:
    base_format = "{:,.0f}"
elif config.base_count_multiplier == 0.000000001:
    base_format = "{:,.2f}"


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
    genstats_headers = {}

    # Init headers for an own separate table
    own_tabl_headers = {}

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


def exist_and_number(data, *metrics):
    return all(isinstance(data.get(m, None), int) or isinstance(data.get(m, None), float) for m in metrics)


def check_duplicate_samples(sample_names, logger):
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
        message = "\n\nDuplicate sample names were found. The last one overwrites previous data." + message
        logger.debug(message)


def order_headers(headers):
    """Produces a shallow copy with ordered single headers."""
    """
    Collected "order_priority" values are stored as groups to resolve
    conflicts, which may arise if several metrics have the same priority.
    The order of headers within the same priority group is not specified.
    """
    indexes = set()
    ordered_headers = {}
    not_ordered_headers = {}
    for metric in headers:
        if "order_priority" in headers[metric] and isinstance(headers[metric]["order_priority"], (int, float)):
            order_priority = headers[metric]["order_priority"]
            indexes.add(order_priority)
            if order_priority not in ordered_headers:
                ordered_headers[order_priority] = []
            ordered_headers[order_priority].append((metric, headers[metric].copy()))
        else:
            not_ordered_headers[metric] = headers[metric].copy()

    # Convert to list and sort in ascending order.
    if indexes:
        indexes = list(indexes)
        indexes.sort()
    else:
        return headers

    output_headers = {}
    for index in indexes:
        for metric in ordered_headers[index]:
            output_headers[metric[0]] = metric[1]

    output_headers.update(not_ordered_headers)
    return output_headers


# STD_TABLE_CONFIGS contains all standard table configurations from:
# https://github.com/MultiQC/MultiQC/blob/main/docs/plots.md#creating-a-table
STD_TABLE_CONFIGS = [
    "namespace",
    "title",
    "description",
    "max",
    "min",
    "ceiling",
    "floor",
    "minrange",
    "scale",
    "bgcols",
    "color",
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
    cleaned_headers = {}
    for metric in headers:
        cleaned_headers[metric] = {
            config: val for config, val in headers[metric].items() if config in STD_TABLE_CONFIGS
        }
    return cleaned_headers


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
    },
    "invalid_file_lines": {
        "coverage_metrics": "\n\nThe lines in files must be:\n"
        "COVERAGE SUMMARY,,<metric>,<value1> or "
        "COVERAGE SUMMARY,,<metric>,<value1>,<value2>\n\n"
        "The following files contain invalid lines:\n",
    },
    "unknown_metrics": {
        "coverage_metrics": "\n\nThe following metrics could not be recognized. "
        "Their headers might be incomplete, hence uninformative and ugly table's columns.\n",
        # ------------- #
        "overall_mean_cov_metrics": "\n\nOnly the 'Average alignment coverage over <file>' metric is supported.\n"
        "The following metrics could not be recognized and are not included in the report.\n",
    },
    "unusual_values": {
        "coverage_metrics": "\n\nAll metrics' values except for int, float and NA are non-standard.\n"
        "The following files contain non-standard values:\n",
        # ------------- #
        "overall_mean_cov_metrics": "\n\nAll metrics' values except for float are non-standard.\n"
        "The following files contain non-standard values:\n",
    },
}


def make_log_report(log_data, logger, module):
    """The only purpose of this function is to create a readable and informative log output
    about found info/warnings/errors, which were found at the time of executing a parser."""

    if "invalid_file_names" in log_data and log_data["invalid_file_names"]:
        if module in DRAGEN_MODULE_TEXTS["invalid_file_names"]:
            log_message = DRAGEN_MODULE_TEXTS["invalid_file_names"][module]
        else:
            log_message = "\n\nThe following files are not valid:\n"

        for root in log_data["invalid_file_names"]:
            log_message += "  " + root + ":\n"
            for file in log_data["invalid_file_names"][root]:
                log_message += "    " + file + "\n"

        logger.warning(log_message + "\n")

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

        logger.debug(log_message + "\n")

    if "unknown_metrics" in log_data and log_data["unknown_metrics"]:
        if module in DRAGEN_MODULE_TEXTS["unknown_metrics"]:
            log_message = DRAGEN_MODULE_TEXTS["unknown_metrics"][module]
        else:
            log_message = "\n\nThe following metrics could not be recognized:\n"

        for metric in log_data["unknown_metrics"]:
            log_message += "  " + metric + "\n"

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

        logger.debug(log_message + "\n")
