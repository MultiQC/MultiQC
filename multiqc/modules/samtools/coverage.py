import copy

from typing import Dict, Union

import logging

from multiqc import BaseMultiqcModule
from multiqc.plots import linegraph, table
from multiqc.plots.table_object import TableConfig

log = logging.getLogger(__name__)


def parse_samtools_coverage(module: BaseMultiqcModule):
    """Find Samtools coverage logs and parse their data"""

    data_by_sample = dict()
    for f in module.find_log_files("samtools/coverage"):
        metrics_by_chrom = parse_single_report(f)
        if len(metrics_by_chrom) > 0:
            if f["s_name"] in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            module.add_data_source(f, section="coverage")
            data_by_sample[f["s_name"]] = metrics_by_chrom

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    # Replace None with actual version if it is available
    module.add_software_version(None)

    # Write parsed report data to a file (restructure first)
    module.write_data_file(data_by_sample, "multiqc_samtools_coverage")

    # Make a table/violin summarising the stats across all regions (mean or sum)
    summary_table(module, data_by_sample)

    # Make a line plot showing coverage stats per region, with a tab switch between stats
    lineplot_per_region(module, data_by_sample)

    # Return the number of logs that were found
    return len(data_by_sample)


def summary_table(module, data_by_sample):
    table_data: Dict[str, Dict[str, float]] = {sname: {} for sname in data_by_sample}
    for sample, d_by_chrom in data_by_sample.items():
        vals = list(d_by_chrom.values())
        table_data[sample]["numreads"] = sum([m["numreads"] for m in vals])
        table_data[sample]["covbases"] = sum([m["covbases"] for m in vals])
        for m in vals:
            m["size"] = m["endpos"] - m["startpos"] + 1
        total_size = sum([m["size"] for m in vals])
        # Average weighted by size. Multiplying by individual weight and dividing by total weight
        table_data[sample]["coverage"] = sum([m["coverage"] * m["size"] for m in vals]) / total_size
        table_data[sample]["meandepth"] = sum([m["meandepth"] * m["size"] for m in vals]) / total_size
        table_data[sample]["meanbaseq"] = sum([m["meanbaseq"] * m["size"] for m in vals]) / total_size
        table_data[sample]["meanmapq"] = sum([m["meanmapq"] * m["size"] for m in vals]) / total_size

    headers: Dict[str, Dict[str, Union[str, bool, int]]] = {
        "numreads": {
            "title": "Reads",
            "description": "Total number of mapped reads",
            "shared_key": "read_count",
            "scale": "RdYlGn",
        },
        "covbases": {
            "title": "Bases",
            "description": "Total number of mapped base pairs",
            "shared_key": "base_count",
            "scale": "Blues",
        },
        "coverage": {
            "title": "Coverage",
            "description": "Percentage of region covered with reads",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "YlGn",
        },
        "meandepth": {
            "title": "Mean depth",
            "description": "Mean depth of coverage",
            "min": 0,
            "suffix": "x",
            "scale": "RdYlGn",
        },
        "meanbaseq": {
            "title": "Mean BQ",
            "description": "Mean base quality",
            "min": 0,
            "scale": "Blues",
        },
        "meanmapq": {
            "title": "Mean MQ",
            "description": "Mean mapping quality",
            "min": 0,
            "max": 60,
            "scale": "RdYlGn",
        },
    }

    module.add_section(
        name="Coverage: global stats",
        anchor="samtools-coverage-table-section",
        description=(
            "Stats parsed from <code>samtools coverage</code> output, and summarized "
            "(added up or weighted-averaged) across all regions."
        ),
        plot=table.plot(
            table_data,
            copy.deepcopy(headers),
            pconfig=TableConfig(
                id="samtools-coverage-table",
                title="Samtools: coverage: Summary",
            ),
        ),
    )

    for h in headers:
        headers[h]["hidden"] = True
    headers["meandepth"]["hidden"] = False
    module.general_stats_addcols(table_data, headers, namespace="coverage")


def lineplot_per_region(module, data_by_sample: Dict):
    tabs: Dict[str, Dict[str, Union[str, int]]] = {
        "numreads": {
            "title": "Mapped reads per region",
            "name": "Reads",
            "ylab": "# reads",
            "tt_suffix": " reads",
            "tt_decimals": 0,
        },
        "covbases": {
            "title": "Mapped base pairs per region",
            "name": "Bases",
            "ylab": "# bases",
            "tt_suffix": " bp",
            "tt_decimals": 0,
        },
        "coverage": {
            "title": "Percentage of region covered with reads",
            "name": "Coverage",
            "ylab": "Coverage %",
            "tt_suffix": "%",
        },
        "meandepth": {
            "title": "Mean depth per region",
            "name": "Mean depth",
            "ylab": "Depth",
            "tt_suffix": "x",
        },
        "meanbaseq": {
            "title": "Mean base quality per region",
            "name": "BQ",
            "ylab": "Base quality",
            "tt_suffix": "",
        },
        "meanmapq": {
            "title": "Mean mapping quality per region",
            "name": "MQ",
            "ylab": "Mapping quality",
            "tt_suffix": "",
        },
    }
    datasets = [
        {
            sample: {chrom: metrics[metric] for chrom, metrics in data_by_sample[sample].items()}
            for sample in data_by_sample
        }
        for metric in tabs
    ]
    data_labels = list(tabs.values())
    for dconfig in data_labels:
        dconfig["title"] = f"Samtools: coverage: {dconfig['title']}"

    module.add_section(
        name="Coverage: stats per region",
        description="Per-region stats parsed from <code>samtools coverage</code> output.",
        anchor="samtools-coverage-section",
        plot=linegraph.plot(
            datasets,
            {
                "id": "samtools-coverage",
                "title": data_labels[0]["title"],
                "xlab": "Region",
                "ylab": data_labels[0]["ylab"],
                "ymin": 0,
                "categories": True,
                "smooth_points": 500,
                "logswitch": True,
                "hide_empty": False,
                "data_labels": data_labels,
            },
        ),
    )


EXPECTED_COLUMNS = [
    "rname",
    "startpos",
    "endpos",
    "numreads",
    "covbases",
    "coverage",
    "meandepth",
    "meanbaseq",
    "meanmapq",
]


def parse_single_report(f) -> Dict[str, Dict[str, Union[int, float]]]:
    """
    Example:
    #rname	startpos	endpos	numreads	covbases	coverage	meandepth	meanbaseq	meanmapq
    oligo_1512_adapters	1	156	10	156	100	7.26282	18.4	29.6
    oligo_741_adapters	1	156	0	0	0	0	0	0

    Returns a dictionary with the contig name (rname) as the key and the rest of the fields as a dictionary
    """
    parsed_data = {}
    lines = f["f"].splitlines()
    expected_header = "#" + "\t".join(EXPECTED_COLUMNS)
    if lines[0] != expected_header:
        logging.warning(f"Expected header for samtools coverage: {expected_header}, got: {lines[0]}")
        return {}

    for idx in range(1, len(lines)):
        line = lines[idx]
        fields = line.strip().split("\t")
        if len(fields) != len(EXPECTED_COLUMNS):
            logging.warning(f"Skipping line with {len(fields)} fields, expected {len(EXPECTED_COLUMNS)}: {line}")
        rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq = fields
        if rname in parsed_data:
            logging.warning(f"Duplicate region found in '{f['s_name']}': {rname}")
            continue
        try:
            parsed_data[rname] = dict(
                startpos=int(startpos),
                endpos=int(endpos),
                numreads=int(numreads),
                covbases=int(covbases),
                coverage=float(coverage),
                meandepth=float(meandepth),
                meanbaseq=float(meanbaseq),
                meanmapq=float(meanmapq),
            )
        except ValueError:
            logging.warning(f"Ignoring invalid line:{idx} for '{f['s_name']}', starting '{line[:6]}'")

    return parsed_data
