"""MultiQC submodule to parse output from GATK BaseRecalibrator"""

import logging
from collections import namedtuple
from itertools import groupby

from multiqc.modules.gatk.utils import parse_report
from multiqc.plots import linegraph, scatter
from multiqc.base_module import BaseMultiqcModule

# Initialise the logger
log = logging.getLogger(__name__)

RecalTableType = namedtuple("RecalTableType", "pre_recalibration post_recalibration".split())
recal_table_type = RecalTableType(0, 1)


def parse_gatk_base_recalibrator(module: BaseMultiqcModule) -> int:
    """Find GATK BaseRecalibrator logs and parse their data"""

    PREFIXES = {
        "gatk": "#:GATKTable:",
        "sentieon": "#:SENTIEON_QCAL_TABLE:",
    }

    REPORT_TABLE_HEADERS = {
        "Arguments:Recalibration argument collection values used in this run": "arguments",
        "Quantized:Quality quantization map": "quality_quantization_map",
        "RecalTable0:": "recal_table_0",
        "RecalTable1:": "recal_table_1",
        "RecalTable2:": "recal_table_2",
    }
    samples_kept = {rt_type: set() for rt_type in recal_table_type}
    data = {
        recal_type: {table_name: {} for table_name in REPORT_TABLE_HEADERS.values()} for recal_type in recal_table_type
    }

    for f in module.find_log_files("gatk/base_recalibrator"):
        # Check that we're not ignoring this sample name
        if module.is_ignore_sample(f["s_name"]):
            continue

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        module.add_software_version(None, f["s_name"])

        # Try GATK first
        parsed_data = parse_report(
            f["f"].splitlines(),
            {PREFIXES["gatk"] + k: v for k, v in REPORT_TABLE_HEADERS.items()},
        )
        if len(parsed_data) == 0:
            parsed_data = parse_report(
                f["f"].splitlines(),
                {PREFIXES["sentieon"] + k: v for k, v in REPORT_TABLE_HEADERS.items()},
            )
        if len(parsed_data) == 0:
            continue

        rt_type = determine_recal_table_type(parsed_data)
        if len(parsed_data) > 0:
            if f["s_name"] in samples_kept[rt_type]:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")
            samples_kept[rt_type].add(f["s_name"])

            module.add_data_source(f, section="base_recalibrator")
            for table_name, sample_tables in parsed_data.items():
                data[rt_type][table_name][f["s_name"]] = sample_tables

    # Filter to strip out ignored sample names. Again.
    for rt_type in recal_table_type:
        for table_name, sample_tables in data[rt_type].items():
            data[rt_type][table_name] = module.ignore_samples(sample_tables)

    n_reports_found = sum([len(samples_kept[rt_type]) for rt_type in recal_table_type])
    if n_reports_found == 0:
        return 0

    log.info(f"Found {n_reports_found} BaseRecalibrator reports")

    # Write data to file
    module.write_data_file(data, "gatk_base_recalibrator")

    add_quality_score_vs_no_of_observations_section(module, data)
    add_reported_vs_empirical_section(module, data)
    return n_reports_found


def add_quality_score_vs_no_of_observations_section(module, data):
    """Add a section for the quality score vs number of observations line plot"""

    sample_data = []
    data_labels = []

    # Loop through the different data types
    for rt_type_name, rt_type in recal_table_type._asdict().items():
        sample_tables = data[rt_type]["quality_quantization_map"]
        if len(sample_tables) == 0:
            continue

        count_data = {}
        pct_data = {}
        for sample, table in sample_tables.items():
            count_data[sample] = {}
            pct_data[sample] = {}

            # Get the total count for this sample
            sample_y_sum = sum(int(y) for y in table["Count"])

            # Collect the data for the plots
            for x, y in zip(table["QualityScore"], table["Count"]):
                # Quality score counts
                count_data[sample][int(x)] = int(y)
                # Quality score percentages
                try:
                    pct_data[sample][int(x)] = float(y) / sample_y_sum
                except ZeroDivisionError:
                    pct_data[sample][int(x)] = 0

        # Append the datasets for this data type
        sample_data.append(count_data)
        sample_data.append(pct_data)

        # Build data label configs for this data type
        data_labels.append({"name": f"{rt_type_name.capitalize().replace('_', '-')} Count", "ylab": "Count"})
        data_labels.append({"name": f"{rt_type_name.capitalize().replace('_', '-')} Percent", "ylab": "Percent"})

    plot = linegraph.plot(
        sample_data,
        pconfig={
            "title": "GATK: Observed Quality Score Counts",
            "id": "gatk-base-recalibrator-quality-scores-plot",
            "xlab": "Observed Quality Score",
            "ylab": "Count",
            "x_decimals": False,
            "data_labels": data_labels,
        },
    )

    # Reported vs empirical quality scores
    module.add_section(
        name="Observed Quality Scores",
        anchor="gatk-base-recalibrator-quality-scores",
        description=(
            "This plot shows the distribution of base quality scores in each sample before and "
            "after base quality score recalibration (BQSR). Applying BQSR should broaden the "
            "distribution of base quality scores."
        ),
        helptext=(
            "For more information see "
            "[the Broad's description of BQSR]"
            "(https://gatkforums.broadinstitute.org/gatk/discussion/44/base-quality-score-recalibration-bqsr)"
            "."
        ),
        plot=plot,
    )


def add_reported_vs_empirical_section(module, data):
    sample_data = []
    data_labels = []

    # Loop through the different data types
    for (
        rt_type_name,
        rt_type,
    ) in recal_table_type._asdict().items():
        # This table appears to be the correct one to use for reported vs empirical
        # https://github.com/broadinstitute/gatk/blob/853b53ec2a3ac2d90d7d82a6c8451e29a34692d2/src/main/resources/org/broadinstitute/hellbender/utils/recalibration/BQSR.R#L148
        sample_tables = data[rt_type]["recal_table_1"]
        if len(sample_tables) == 0:
            continue

        reported_empirical = {}
        for sample, table in sample_tables.items():
            reported_empirical[sample] = []
            table_rows = [dict(zip(table, r)) for r in zip(*table.values())]
            table_rows.sort(key=lambda r: r["QualityScore"])
            for reported, group in groupby(table_rows, lambda r: r["QualityScore"]):
                g = list(group)
                reported_empirical[sample].append(
                    {
                        "x": int(reported),
                        "y": sum(float(r["EmpiricalQuality"]) for r in g) / len(g) if len(g) > 0 else 0,
                    }
                )

        sample_data.append(reported_empirical)

        # Build data label configs for this data type
        data_labels.append(
            {"name": f"{rt_type_name} Reported vs. Empirical Quality", "ylab": "Empirical quality score"}
        )

    plot = scatter.plot(
        sample_data,
        pconfig={
            "title": "Reported vs. Empirical Quality",
            "id": "gatk-base-recalibrator-reported-empirical-plot",
            "xlab": "Reported quality score",
            "ylab": "Empirical quality score",
            "x_decimals": False,
            "data_labels": data_labels,
            "xmin": 0,
            "ymin": 0,
            "square": True,
        },
    )

    module.add_section(
        name="Reported Quality vs. Empirical Quality",
        anchor="gatk-base-recalibrator-reported-empirical",
        description="Plot shows the reported quality score vs the empirical quality score.",
        plot=plot,
    )


def determine_recal_table_type(parsed_data):
    arguments = dict(zip(parsed_data["arguments"]["Argument"], parsed_data["arguments"]["Value"]))
    if arguments["recalibration_report"] == "null":
        return recal_table_type.pre_recalibration
    else:
        return recal_table_type.post_recalibration
