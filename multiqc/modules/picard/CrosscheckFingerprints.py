""" MultiQC submodule to parse output from Picard CrosscheckFingerprints """

import logging
import re

from collections import OrderedDict
from csv import DictReader
from distutils.util import strtobool
from itertools import chain
from itertools import groupby
from multiqc import config
from multiqc.plots import table

# Initialize the logger
log = logging.getLogger(__name__)

# This is a subset, the rest of the fields are self descriptive
FIELD_DESCRIPTIONS = {
    "LEFT_GROUP_VALUE": "The name of the left data-type group.",
    "RIGHT_GROUP_VALUE": "The name of the right data-type grup.",
    "RESULT": "The categorical result of comparing the calculated LOD score against the threshold.",
    "DATA_TYPE": "The datatype used for the comparison.",
    "LOD_SCORE": "Log10 of the probability that the samples come from the same individual.",
    "LOD_SCORE_TUMOR_NORMAL": "LOD score with the assumption that Left is a Tumor.",
    "LOD_SCORE_NORMAL_TUMOR": "LOD score with the assumption that Right is a Tumor.",
    "LOD_THRESHOLD": "The LOD threshold used for this pairwise comparison.",
    "TUMOR_AWARENESS": "Whether or not this pairwise comparison was flagged for tumor awareness",
}


def parse_reports(self):
    """Find Picard CrosscheckFingerprints reports and parse their data.

    Stores the data in "Sample/Group - Sample/Group" groups since CrosscheckFingerprints
    does pairwise comparisons between samples at the level selected by `--CROSSCHECK_BY`.
    """

    self.picard_CrosscheckFingerprints_data = dict()

    # Go through logs and find Metrics
    for f in self.find_log_files("picard/crosscheckfingerprints", filehandles=True):
        # Parse an individual CrosscheckFingerprints Report
        (metrics, comments) = _take_till(f["f"], lambda line: line.startswith("#") or line == "\n")
        header = next(metrics).rstrip("\n").split("\t")
        if not "LEFT_GROUP_VALUE" in header:
            # Not a CrosscheckFingerprints Report
            continue
        reader = DictReader(metrics, fieldnames=header, delimiter="\t")
        # Parse out the tumor awareness option and the lod threshold setting if possible
        (tumor_awareness, lod_threshold) = _parse_cli(comments[1])
        for row in reader:
            # Decide what to use as the `sample_name` for this row
            name = _get_name(row)

            # Set the cli options of interest for this file
            row["LOD_THRESHOLD"] = lod_threshold
            row["TUMOR_AWARENESS"] = tumor_awareness
            self.picard_CrosscheckFingerprints_data[name] = row

    # For each sample, flag if any comparisons don't start with "Expected"
    # A sample that does not have all "Expected" will show as `False` and be Red
    general_stats_data = _create_general_stats_data(self.picard_CrosscheckFingerprints_data)
    self.general_stats_addcols(
        general_stats_data,
        {
            "Crosschecks All Expected": {
                "title": "Crosschecks All Expected",
                "description": "All results for samples CrosscheckFingerprints were as expected.",
                "scale": "RdYlGn",
            }
        },
    )

    # Add a table section to the report
    self.add_section(
        name="Crosscheck Fingerprints",
        anchor="picard-crosscheckfingerprints",
        description="Pairwise identity checking betwen samples and groups.",
        helptext="""
        Checks that all data in the set of input files comes from the same individual, based on the selected group granularity.

        """,
        plot=table.plot(
            self.picard_CrosscheckFingerprints_data,
            _get_table_headers(self.picard_CrosscheckFingerprints_data),
        ),
    )

    return len(self.picard_CrosscheckFingerprints_data)


def _take_till(iterator, fn):
    """Take from an iterator till `fn` returns false.

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

    return (chain([val], iterator), headers)


def _parse_cli(line):
    """Parse the Picard CLI invocation that is stored in the header section of the file."""
    tumor_awareness_regex = r"CALCULATE_TUMOR_AWARE_RESULTS=(\w+)"
    lod_threshold_regex = r"LOD_THRESHOLD=(\S+)"

    tumor_awareness = None
    lod_threshold = None

    tumor_awareness_match = re.search(tumor_awareness_regex, line)
    if tumor_awareness_match is not None:
        tumor_awareness = strtobool(tumor_awareness_match.group(1))

    lod_threshold_match = re.search(lod_threshold_regex, line)
    if lod_threshold_match is not None:
        lod_threshold = float(lod_threshold_match.group(1))

    return (tumor_awareness, lod_threshold)


def _get_name(row):
    """ Make a name for the row based on the Samples and Groups present """
    if (
        row["LEFT_GROUP_VALUE"] == row["LEFT_SAMPLE"]
        and row["RIGHT_GROUP_VALUE"] == row["RIGHT_SAMPLE"]
    ):
        return "{} - {}".format(row["LEFT_SAMPLE"], row["RIGHT_SAMPLE"])
    else:
        return "{}/{} - {}/{}".format(
            row["LEFT_SAMPLE"],
            row["LEFT_GROUP_VALUE"],
            row["RIGHT_SAMPLE"],
            row["RIGHT_GROUP_VALUE"],
        )


def _get_table_headers(data):
    """ Create the headers config """

    picard_config = getattr(config, "picard_config", {})
    crosscheckfingerprints_table_cols = picard_config.get("CrosscheckFingerprints_table_cols")
    crosscheckfingerprints_table_cols_hidden = picard_config.get(
        "CrosscheckFingerprints_table_cols_hidden"
    )

    if not crosscheckfingerprints_table_cols:
        crosscheckfingerprints_table_cols = [
            "RESULT",
            "DATA_TYPE",
            "LOD_THRESHOLD",
            "LOD_SCORE",
        ]
    if not crosscheckfingerprints_table_cols_hidden:
        crosscheckfingerprints_table_cols_hidden = [
            "LEFT_GROUP_VALUE",
            "RIGHT_GROUP_VALUE",
            "LEFT_RUN_BARCODE",
            "LEFT_LANE",
            "LEFT_MOLECULAR_BARCODE_SEQUENCE",
            "LEFT_LIBRARY",
            "LEFT_SAMPLE",
            "LEFT_FILE",
            "RIGHT_RUN_BARCODE",
            "RIGHT_LANE",
            "RIGHT_MOLECULAR_BARCODE_SEQUENCE",
            "RIGHT_LIBRARY",
            "RIGHT_SAMPLE",
            "RIGHT_FILE",
        ]

        # Add the Tumor/Normal LOD scores if any pair had the tumor_awareness flag set
        if any(row["TUMOR_AWARENESS"] for row in data.values()):
            crosscheckfingerprints_table_cols += [
                "LOD_SCORE_TUMOR_NORMAL",
                "LOD_SCORE_NORMAL_TUMOR",
            ]
        else:
            crosscheckfingerprints_table_cols_hidden += [
                "LOD_SCORE_TUMOR_NORMAL",
                "LOD_SCORE_NORMAL_TUMOR",
            ]

        headers = OrderedDict()
        for h in FIELD_DESCRIPTIONS:
            # Skip anything not set to visible
            if h not in crosscheckfingerprints_table_cols:
                continue

            # Set up the configuration for the column
            h_title = h.replace("_", " ").strip().lower().capitalize()
            headers[h] = {
                "title": h_title,
                "description": FIELD_DESCRIPTIONS.get(h),
                "namespace": "CrosscheckFingerprints",
            }

            # Add appropriate colors for LOD scores
            if h.startswith("LOD"):
                headers[h]["scale"] = "RdYlGn"

            if h in crosscheckfingerprints_table_cols_hidden:
                headers[h]["hidden"] = True

        return headers


def _create_general_stats_data(in_data):
    """Look at the LEFT_SAMPLE fields and determine if there are any pairs for that samples
    that don't have a RESULT that startswith EXPECTED.
    """
    out_data = dict()
    flattened = (row for row in in_data.values())
    sorted_by_left_sample = sorted(flattened, key=lambda r: r["LEFT_SAMPLE"])

    for group, values in groupby(sorted_by_left_sample, key=lambda r: r["LEFT_SAMPLE"]):
        out_data[group] = {
            # NB: Must coerce bool to str or else plot.table turns it into a float
            "Crosschecks All Expected": str(
                all(v["RESULT"].startswith("EXPECTED") for v in values)
            )
        }

    return out_data
