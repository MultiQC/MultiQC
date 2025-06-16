import csv
import io
import logging
from typing import Dict, Union, cast, Any, List
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table
from multiqc.plots.table_object import TableConfig, ColumnMeta, ColumnDict
from pprint import pprint

# No DOI available – the data parsed here is custom and not from a published source
# self.add_software_version not needed – no external software version applies

log = logging.getLogger(__name__)


def parse_csv_metrics(fh, filename: str) -> Dict[str, Dict[str, Union[float, int]]]:
    data: Dict[str, Dict[str, Union[int, float]]] = {}
    sample_name = filename.replace(".csv", "")

    if isinstance(fh, io.BytesIO):
        fh = io.TextIOWrapper(fh, encoding="latin-1")
    else:
        try:
            fh = io.TextIOWrapper(fh.buffer, encoding="latin-1")
        except AttributeError:
            pass

    reader = csv.DictReader(fh)
    for row in reader:
        metric_name = row.get("Metric Name", "").strip()
        value = row.get("Metric Value", "").strip()
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
    data: Dict[str, Dict[str, Union[int, float]]], headers: Dict[str, ColumnDict]
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
        # self.doi = None  # No DOI – this module parses custom CSV data

        data_by_sample: Dict[str, Dict[str, Union[int, float]]] = {}

        for f in self.find_log_files("custom_csv", filehandles=True):
            d_name = f["root"]
            # Clean the sample name for MultiQC standards
            sample_name = self.clean_s_name(d_name, f)

            self.add_software_version(None, sample_name)

            # Parse the file using sample_name
            parsed = parse_csv_metrics(f["f"], sample_name)
            if not parsed:
                continue

            for sample_id, sample_metrics in parsed.items():
                data_by_sample[sample_id] = sample_metrics

            self.add_data_source(f)

        if not data_by_sample:
            raise ModuleNoSamplesFound

        self.write_data_file(data_by_sample, "multiqc_csv")
        color_palette: List[str] = [
            "OrRd",
            "BuPu",
            "Oranges",
            "BuGn",
            "PuBu",
            "PuRd",
            "PuBuGn",
            "Reds",
            "RdPu",
            "Greens",
            "Greys",
            "Purples",
            "GnBu",
        ]

        headers: Dict[str, ColumnDict] = {}
        color_index = 0
        for sample_metrics in data_by_sample.values():
            for metric in sample_metrics:
                if metric not in headers:
                    headers[metric] = {
                        "title": metric,
                        "description": f"Parsed metric: {metric}",
                        "scale": color_palette[color_index % len(color_palette)],
                    }
                    color_index += 1

        # Plot the custom summary table
        self.add_section(
            name="CSV Summary",
            anchor="csv_summary",
            description="Metrics parsed from the input CSV file.",
            plot=table.plot(
                data_by_sample,
                headers,
                pconfig=TableConfig(id="csv_summary_table", title="CSV Summary Metrics"),
            ),
        )

        # Configure which general stats columns are hidden or visible
        general_stats_headers = {k: dict(v) for k, v in headers.items()}
        for col in general_stats_headers.values():
            col["hidden"] = True
        general_stats_headers["Cells"]["hidden"] = False
        general_stats_headers["Number of reads"]["hidden"] = False
        self.general_stats_addcols(format_general_stats(data_by_sample, headers), general_stats_headers)
