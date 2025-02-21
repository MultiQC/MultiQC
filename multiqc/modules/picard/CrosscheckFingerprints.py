"""MultiQC submodule to parse output from Picard CrosscheckFingerprints"""

import logging
import re
from collections import OrderedDict, defaultdict
from csv import DictReader
from itertools import chain, groupby

from multiqc import config
from multiqc.plots import table, heatmap
from multiqc.utils.util_functions import strtobool

# Initialize the logger
log = logging.getLogger(__name__)

# This is a subset, the rest of the fields are self-descriptive
FIELD_DESCRIPTIONS = {
    "LEFT_SAMPLE": "The name of the left sample.",
    "LEFT_GROUP_VALUE": "The name of the left data-type group.",
    "RIGHT_SAMPLE": "The name of the right sample.",
    "RIGHT_GROUP_VALUE": "The name of the right data-type group.",
    "RESULT": "The categorical result of comparing the calculated LOD score against the threshold.",
    "DATA_TYPE": "The datatype used for the comparison.",
    "LOD_SCORE": "Log10 of the probability that the samples come from the same individual.",
    "LOD_SCORE_TUMOR_NORMAL": "LOD score with the assumption that Left is a Tumor.",
    "LOD_SCORE_NORMAL_TUMOR": "LOD score with the assumption that Right is a Tumor.",
    "LOD_THRESHOLD": "The LOD threshold used for this pairwise comparison.",
    "TUMOR_AWARENESS": "Whether this pairwise comparison was flagged for tumor awareness",
}


def parse_reports(module):
    """
    Find Picard CrosscheckFingerprints reports and parse their data.

    Stores the data in "Sample/Group - Sample/Group" groups since CrosscheckFingerprints
    does pairwise comparisons between samples at the level selected by `--CROSSCHECK_BY`.
    """

    row_by_number = dict()
    found_reports = []

    # Go through logs and find Metrics
    row_number = 0
    for f in module.find_log_files("picard/crosscheckfingerprints", filehandles=True):
        # Parse an individual CrosscheckFingerprints Report
        (metrics, comments) = _take_till(f["f"], lambda line: line.startswith("#") or line == "\n")
        header = next(metrics).rstrip("\n").split("\t")
        if "LEFT_GROUP_VALUE" not in header:
            # Not a CrosscheckFingerprints Report
            continue
        # Parse out the tumor awareness option and the lod threshold setting if possible
        tumor_awareness, lod_threshold = _parse_cli(comments[1])
        reader: DictReader = DictReader(metrics, fieldnames=header, delimiter="\t")
        lines = []
        for row in reader:
            if not module.is_ignore_sample(row["LEFT_SAMPLE"]) and not module.is_ignore_sample(row["RIGHT_SAMPLE"]):
                lines.append((row_number, row))
            row_number += 1
        if not lines:
            continue
        found_reports.append(f["s_name"])
        for row_number, row in lines:
            # Clean the sample names
            row["LEFT_SAMPLE"] = module.clean_s_name(row["LEFT_SAMPLE"], f)
            row["LEFT_GROUP_VALUE"] = module.clean_s_name(row["LEFT_GROUP_VALUE"], f)
            row["RIGHT_SAMPLE"] = module.clean_s_name(row["RIGHT_SAMPLE"], f)
            row["RIGHT_GROUP_VALUE"] = module.clean_s_name(row["RIGHT_GROUP_VALUE"], f)

            row["RESULT"] = row["RESULT"].capitalize().replace("_", " ")

            # Set the cli options of interest for this file
            row["LOD_THRESHOLD"] = lod_threshold
            row["TUMOR_AWARENESS"] = tumor_awareness
            row_by_number[row_number] = row

            try:
                row["LOD_SCORE"] = float(row["LOD_SCORE"])
            except ValueError:
                row["LOD_SCORE"] = None

            # Add BEST_MATCH and BEST_MATCH_LOD if the result is unexpected or inconclusive
            if row["RESULT"].startswith("Unexpected") or row["RESULT"] == "Inconclusive":
                # Find the best match for the left sample
                left_sample = row["LEFT_SAMPLE"]
                best_match = max(
                    (r for j, r in lines if r["LEFT_SAMPLE"] == left_sample), key=lambda r: float(r["LOD_SCORE"])
                )
                row["BEST_MATCH"] = best_match["RIGHT_SAMPLE"]
                row["BEST_MATCH_LOD"] = float(best_match["LOD_SCORE"])

            module.add_data_source(f, section="CrosscheckFingerprints")

    # Only add sections if we found data
    if not found_reports:
        return set()

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write data to file
    module.write_data_file(row_by_number, f"{module.anchor}_crosscheckfingerprints")

    # Add a per-sample table
    status_by_sample = dict()
    sorted_by_left_sample = sorted(row_by_number.values(), key=lambda r: r["LEFT_SAMPLE"])
    for left_sample, values in groupby(sorted_by_left_sample, key=lambda r: r["LEFT_SAMPLE"]):
        values = list(values)
        if values:
            status = "All expected"
            if all(v["RESULT"].startswith("Unexpected") for v in values):
                status = "All unexpected"
            elif any(v["RESULT"].startswith("Unexpected") for v in values):
                status = "Some unexpected"
            elif any(v["RESULT"] == "Inconclusive" for v in values):
                status = "Some inconclusive"
            status_by_sample[left_sample] = {"Crosschecks": status}
            best_match = max(values, key=lambda v: v["LOD_SCORE"])
            status_by_sample[left_sample]["Best match"] = best_match["RIGHT_SAMPLE"]
            status_by_sample[left_sample]["Best match LOD"] = best_match["LOD_SCORE"]

    sample_table_headers = {
        "Crosschecks": {
            "title": "Crosschecks",
            "description": "Flags if there are any unexpected or inconclusive matches.",
            "cond_formatting_rules": {
                "pass": [{"s_eq": "All expected"}],
                "warn": [{"s_eq": "Some inconclusive"}],
                "fail": [{"s_eq": "Unexpected"}, {"s_eq": "Some unexpected"}],
            },
        },
        "Best match": {
            "title": "Best match",
            "description": "The sample name with the highest LOD for this sample.",
        },
        "Best match LOD": {
            "title": "Best match LOD",
            "description": "The LOD score for the best matching sample.",
        },
    }

    def _order_status(status: str) -> int:
        # sorting order to show the bad results first
        if "all unexpected" in status.lower():
            return 0
        if "unexpected" in status.lower():
            return 1
        if "inconclusive" in status.lower():
            return 2
        return 3

    module.add_section(
        name="Crosscheck Fingerprints",
        anchor=f"{module.anchor}-crosscheckfingerprints-sample-table-section",
        description="Summary for each sample, showing the overall status along with the best matching sample name.",
        plot=table.plot(
            dict(sorted(status_by_sample.items(), key=lambda kv: _order_status(kv[1]["Crosschecks"]))),
            sample_table_headers,
            pconfig={
                "namespace": module.name,
                "id": f"{module.anchor}-crosscheckfingerprints-sample-table",
                "title": f"{module.name}: Crosscheck Fingerprints: Samples",
                "no_violin": True,
                "sort_rows": False,
            },
        ),
    )
    for h in sample_table_headers.values():
        h["hidden"] = True
    module.general_stats_addcols(
        status_by_sample,
        sample_table_headers,
        namespace="CrosscheckFingerprints",
    )

    # Heatmap of the LOD scores for each pairwise comparison
    heatmap_data = defaultdict(dict)
    for row in row_by_number.values():
        left = row["LEFT_SAMPLE"]
        right = row["RIGHT_SAMPLE"]
        heatmap_data[left][right] = row["LOD_SCORE"]

    module.add_section(
        name="Crosscheck Fingerprints: Heatmap",
        anchor=f"{module.anchor}-crosscheckfingerprints-heatmap-section",
        description="Pairwise identity checking between samples and groups: heatmap of LOD scores.",
        plot=heatmap.plot(
            heatmap_data,
            xcats=list(heatmap_data.keys()),
            pconfig={
                "id": f"{module.anchor}-crosscheckfingerprints-heatmap",
                "title": f"{module.name}: Crosscheck Fingerprints",
                "zlab": "LOD score",
                "xcats_samples": True,
                "ycats_samples": True,
                "square": True,
                "legend": False,
                "reverse_colors": True,
            },
        ),
    )

    # If the number of rows is > 100, do not show pairs with "Expected" status
    warning = ""
    if len(row_by_number) > 100:
        warning = (
            f"Note that there are too many pairwise comparisons to show in table "
            f" ({len(row_by_number)} > 100), so only unexpected or inconclusive pairs are shown."
        )
        row_by_number = {k: v for k, v in row_by_number.items() if not v["RESULT"].startswith("Expected")}

    module.add_section(
        name="Crosscheck Fingerprints: Pairwise Table",
        anchor=f"{module.anchor}-crosscheckfingerprints-table-section",
        description="Pairwise identity checking between samples and groups." + (f"<br>{warning}" if warning else ""),
        helptext="""
        Checks that all data in the set of input files comes from the same individual, based on the selected group granularity.
        """,
        plot=table.plot(
            data=dict(sorted(row_by_number.items(), key=lambda kv: (_order_status(kv[1]["RESULT"]), kv[0]))),
            headers=_get_table_headers(row_by_number),
            pconfig={
                "namespace": module.name,
                "id": f"{module.anchor}-crosscheckfingerprints-table",
                "title": f"{module.name}: Crosscheck Fingerprints",
                "save_file": True,
                "col1_header": "ID",
                "no_violin": True,
                "sort_rows": False,
            },
        ),
    )

    return found_reports


def _take_till(iterator, fn):
    """
    Take from an iterator till `fn` returns false.

    Returns the iterator with the value that caused false at the front, and all the lines skipped till then as a list.
    """
    headers = []
    try:
        val = next(iterator)
        while fn(val):
            headers.append(val)
            val = next(iterator)
    except StopIteration:
        return ()

    return chain([val], iterator), headers


def _parse_cli(line):
    """Parse the Picard CLI invocation that is stored in the header section of the file."""
    tumor_awareness_regex = r"CALCULATE_TUMOR_AWARE_RESULTS(\s|=)(\w+)"
    lod_threshold_regex = r"LOD_THRESHOLD(\s|=)(\S+)"

    tumor_awareness = None
    lod_threshold = None

    tumor_awareness_match = re.search(tumor_awareness_regex, line)
    if tumor_awareness_match is not None:
        tumor_awareness = strtobool(tumor_awareness_match.group(2))

    lod_threshold_match = re.search(lod_threshold_regex, line)
    if lod_threshold_match is not None:
        lod_threshold = float(lod_threshold_match.group(2))

    return tumor_awareness, lod_threshold


def _get_table_headers(data_by_sample):
    """Create the headers config"""

    table_cols = [
        "RESULT",
        "DATA_TYPE",
        "LOD_THRESHOLD",
        "LOD_SCORE",
    ]
    table_cols_hidden = [
        "LEFT_RUN_BARCODE",
        "LEFT_LANE",
        "LEFT_MOLECULAR_BARCODE_SEQUENCE",
        "LEFT_LIBRARY",
        "LEFT_FILE",
        "RIGHT_RUN_BARCODE",
        "RIGHT_LANE",
        "RIGHT_MOLECULAR_BARCODE_SEQUENCE",
        "RIGHT_LIBRARY",
        "RIGHT_FILE",
        "DATA_TYPE",
        "LOD_THRESHOLD",
    ]

    # Allow customisation from the MultiQC config
    picard_config = getattr(config, "picard_config", {})
    table_cols = picard_config.get("CrosscheckFingerprints_table_cols", table_cols)
    table_cols_hidden = picard_config.get("CrosscheckFingerprints_table_cols_hidden", table_cols_hidden)

    # Add the Tumor/Normal LOD scores if any pair had the tumor_awareness flag set
    if any(row["TUMOR_AWARENESS"] for row in data_by_sample.values()):
        table_cols += [
            "LOD_SCORE_TUMOR_NORMAL",
            "LOD_SCORE_NORMAL_TUMOR",
        ]
    else:
        table_cols_hidden += [
            "LOD_SCORE_TUMOR_NORMAL",
            "LOD_SCORE_NORMAL_TUMOR",
        ]

    # Add Left and Right Sample names / groups, keeping it as minimal as possible
    def sample_group_are_same(x):
        return x["LEFT_SAMPLE"] == x["LEFT_GROUP_VALUE"] and x["RIGHT_SAMPLE"] == x["RIGHT_GROUP_VALUE"]

    if all(sample_group_are_same(values) for values in data_by_sample.values()):
        table_cols = [
            "LEFT_SAMPLE",
            "RIGHT_SAMPLE",
        ] + table_cols
        table_cols_hidden += ["LEFT_GROUP_VALUE", "RIGHT_GROUP_VALUE"]
    else:
        table_cols = [
            "LEFT_SAMPLE",
            "LEFT_GROUP_VALUE",
            "RIGHT_SAMPLE",
            "RIGHT_GROUP_VALUE",
        ] + table_cols

    headers = OrderedDict()
    for h in FIELD_DESCRIPTIONS:
        # Skip anything not set to visible
        if h not in table_cols:
            continue

        # Set up the configuration for the column
        h_title = (
            h.replace("_", " ")
            .strip()
            .lower()
            .capitalize()
            .replace("Lod score", "LOD")
            .replace("tumor normal", "T/N")
            .replace("normal tumor", "N/T")
        )
        headers[h] = {
            "title": h_title,
            "description": FIELD_DESCRIPTIONS.get(h),
            "namespace": "CrosscheckFingerprints",
            "scale": False,
        }

        # Rename Result to be a longer string so the table formats more nicely
        if h == "RESULT":
            headers[h]["title"] = "Match"
            headers[h]["cond_formatting_rules"] = {
                "pass": [{"s_contains": "Expected"}],
                "warn": [{"s_eq": "Inconclusive"}],
                "fail": [{"s_contains": "Unexpected"}],
            }

        # Add appropriate colors for LOD scores
        if h.startswith("LOD"):
            headers[h]["scale"] = "RdYlGn"
            headers[h]["shared_key"] = "LOD"
            headers[h]["bars_zero_centrepoint"] = True

        if h in table_cols_hidden:
            headers[h]["hidden"] = True

    # Add a new column showing the best match for the left sample if there are unexpected samples
    if "RESULT" in table_cols and "LEFT_SAMPLE" in table_cols and "RIGHT_SAMPLE" in table_cols:
        headers["BEST_MATCH"] = {
            "title": "Best match",
            "description": "The sample name with the highest LOD for the left sample.",
            "namespace": "CrosscheckFingerprints",
            "hidden": False,
        }
        headers["BEST_MATCH_LOD"] = {
            "title": "Best match LOD",
            "description": "The LOD score for the best matching sample for the left sample.",
            "namespace": "CrosscheckFingerprints",
            "hidden": False,
        }

    return headers
