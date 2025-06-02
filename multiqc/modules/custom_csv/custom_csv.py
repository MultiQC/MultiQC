import csv
import logging
from typing import Dict, Union, cast, Any
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import TableConfig, ColumnMeta, ColumnDict
from pprint import pprint

log = logging.getLogger(__name__)


def parse_csv_metrics(fh, filename: str) -> Dict[str, Dict[str, Union[float, int]]]:
    data: Dict[str, Dict[str, Union[int, float]]] = {}
    sample_name = filename.replace(".csv", "")
    reader = csv.DictReader(fh)
    for row in reader:
        metric_name = row.get("Metric Name", "").strip()
        value = row.get("Metric Value", "").strip()
        # s_name = row.get("Group Name", "").strip() or "global"
        if not metric_name or not value:
            continue
        value = value.replace(",", "").replace("%", "")
        try:
            val = float(value)
        except ValueError:
            continue
        if sample_name not in data:
            data[sample_name] = {}
        data[sample_name][metric_name] = val
    return data


def format_general_stats( 
    data: Dict[str, Dict[str, Union[int, float]]],
    headers: Dict[str, ColumnDict]
) -> Dict[str, Dict[str, Union[int, float, str, bool]]]:
    formatted: Dict[str, Dict[str, Union[int, float, str, bool]]] = {}
    for sample, metrics in data.items():
        formatted[sample] = {}
        for k, v in metrics.items():
            formatted[sample][k] = v
    return formatted




class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="CSVFlex",
            anchor="custom_csv",
            info="Parses a CSV file to extract summary metrics and display them in MultiQC.",
        )
        log.info("Hello World!")

        data_by_sample: Dict[str, Dict[str, Union[int, float]]] = {}

        for f in self.find_log_files("custom_csv", filehandles=True):
            parsed = parse_csv_metrics(f["f"], f["fn"])
            pprint(parsed)
            log.info(parsed)
            if not parsed:
                continue
            data_by_sample = parsed
            self.add_data_source(f)

        if not data_by_sample:
            raise ModuleNoSamplesFound

        self.write_data_file(data_by_sample, "multiqc_csv")
        log.info(data_by_sample)

        # headers: Dict[str, ColumnDict] = {
        #     "Cells": {
        #         "title": "Cell Count",
        #         "scale": "Greens",
        #         "description": "Number of cells identified in the sample",
        #         # "rid": "cells",
        #         # "clean_rid": "cells"
        #     },
        #     "Confidently mapped reads in cells": {
        #         "title": "% Reads in Cells",
        #         "suffix": "%",
        #         "description": "Percentage of confidently mapped reads that are within cells",
        #         # "rid": "mapped_reads_in_cells",
        #         # "clean_rid": "mapped_reads_in_cells"
        #     },
        #     "Mean reads per cell": {
        #         "title": "Reads/Cell",
        #         "scale": "Blues",
        #         "description": "Average number of reads per identified cell",
        #         # "rid": "mean_reads_per_cell",
        #         # "clean_rid": "mean_reads_per_cell"
        #     },
        #     "Number of reads": {
        #         "title": "# Reads",
        #         "scale": "Oranges",
        #         "description": "Total number of sequencing reads",
        #         # "rid": "number_of_reads",
        #         # "clean_rid": "number_of_reads"
        #     },
        #     "Valid barcodes": {
        #         "title": "Valid Barcodes",
        #         "suffix": "%",
        #         "description": "Percentage of barcodes deemed valid",
        #         # "rid": "valid_barcodes",
        #         # "clean_rid": "valid_barcodes"
        #     },
        #     "Q30 RNA read": {
        #         "title": "Q30 RNA",
        #         "suffix": "%",
        #         "description": "Percentage of RNA bases with a quality score above Q30",
        #         # "rid": "q30_rna",
        #         # "clean_rid": "q30_rna"
        #     }
        # }

        headers: Dict[str, ColumnDict] = {}
        for sample_metrics in data_by_sample.values():
            for metric in sample_metrics:
                if metric not in headers:
                    headers[metric] = {
                        "title": metric,
                        "description": f"Parsed metric: {metric}"
                    }



        # Plot the custom summary table
        self.add_section(
            name="CSV Summary",
            anchor="csv_summary",
            description="Metrics parsed from the input CSV file.",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig=TableConfig(
                    id="csv_summary_table",
                    title="CSV Summary Metrics"
                ),
            ),
        )

        # Configure which general stats columns are hidden or visible
        general_stats_headers = {k: dict(v) for k, v in headers.items()}
        for col in general_stats_headers.values():
            col["hidden"] = True
        general_stats_headers["Cells"]["hidden"] = False
        general_stats_headers["Number of reads"]["hidden"] = False

        # Format stats for general stats table
        # log.info(data_by_sample)
        formatted_stats = format_general_stats(data_by_sample, headers)
        # pprint(format_general_stats((data_by_sample, headers)))

        print("FINAL general_stats DATA:")
        pprint(formatted_stats)

        self.general_stats_addcols(
            format_general_stats(data_by_sample, headers),
            general_stats_headers
        )
