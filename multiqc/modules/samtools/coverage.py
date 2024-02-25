import copy

from typing import Dict, Union

import logging

from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


class CoverageReportMixin:
    def parse_samtools_coverage(self):
        """Find Samtools coverage logs and parse their data"""

        data_by_sample = dict()
        for f in self.find_log_files("samtools/coverage"):
            metrics_by_chrom = parse_single_report(f)
            if len(metrics_by_chrom) > 0:
                if f["s_name"] in data_by_sample:
                    log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="coverage")
                data_by_sample[f["s_name"]] = metrics_by_chrom

        # Filter to strip out ignored sample names
        data_by_sample = self.ignore_samples(data_by_sample)
        if len(data_by_sample) == 0:
            return 0

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Write parsed report data to a file (restructure first)
        self.write_data_file(data_by_sample, "multiqc_samtools_coverage")

        # Make a table/violin summarising the stats across all regions (mean or sum)
        self.summary_table(data_by_sample)

        # Make a line plot showing coverage stats per region, with a tab switch between stats
        self.lineplot_per_region(data_by_sample)

        # Return the number of logs that were found
        return len(data_by_sample)

    def summary_table(self, data_by_sample):
        table_data = {sname: {} for sname in data_by_sample}
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

        headers = {
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

        self.add_section(
            name="Coverage: global stats",
            anchor="samtools-coverage-table-section",
            description=(
                "Stats parsed from <code>samtools coverage</code> output, and summarized "
                "(added up or weighted-averaged) across all regions."
            ),
            plot=table.plot(
                table_data,
                copy.deepcopy(headers),
                pconfig={
                    "id": "samtools-coverage-table",
                    "title": "Samtools Coverage: Summary",
                },
            ),
        )

        for h in headers:
            headers[h]["hidden"] = True
        headers["meandepth"]["hidden"] = False
        self.general_stats_addcols(table_data, headers, namespace="coverage")

    def lineplot_per_region(self, data_by_sample):
        tabs = {
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
            dconfig["title"] = "Samtools coverage: " + dconfig["title"]

        self.add_section(
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
                    "hide_zero_cats": False,
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

    for line in lines[1:]:
        fields = line.strip().split("\t")
        if len(fields) != len(EXPECTED_COLUMNS):
            logging.warning(f"Skipping line with {len(fields)} fields, expected {len(EXPECTED_COLUMNS)}: {line}")
        rname, startpos, endpos, numreads, covbases, coverage, meandepth, meanbaseq, meanmapq = fields
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

    return parsed_data
