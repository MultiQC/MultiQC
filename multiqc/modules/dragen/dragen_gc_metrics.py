import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


class DragenGcMetrics(BaseMultiqcModule):
    """Not to be confused with DragenFastqcGcMetrics"""

    NAMESPACE = "GC Metrics"

    def add_gc_metrics_hist(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/gc_metrics"):
            data = parse_gc_metrics_file(f)
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="gc_metrics")
            data_by_sample[s_name] = data

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        # Only plot data, don't want to write this to a file
        # (can do so with --export-plots already)
        # self.write_data_file(data_by_sample, "dragen_gc_metrics")

        hist_data = DragenGcMetrics.__get_normalized_gc_data(data_by_sample)
        smooth_points = 300
        self.add_section(
            name="GC Bias Histogram",
            anchor="dragen-gc-bias-hist",
            description="""
                A histogram of the normalized coverage vs GC content.  This shows how GC
                content in the genome impacts sequencing coverage.
                """,
            plot=linegraph.plot(
                hist_data,
                {
                    "id": "gc-bias-hist",
                    "title": "Dragen: GC Bias Histogram",
                    "ylab": "Normalized Coverage",
                    "xlab": "% GC",
                    "ymin": 0,
                    "xmin": 0,
                    "tt_label": "<b>{point.x} % GC</b>: {point.y} Normalized coverage",
                    "smooth_points": smooth_points,
                    "namespace": DragenGcMetrics.NAMESPACE,
                },
            ),
        )

        table_data = DragenGcMetrics.__get_summary_gc_data(data_by_sample)
        self.add_section(
            name="GC Metrics Summary",
            anchor="dragen-gc-metrics-summary",
            description="""
            Summary GC metrics shown on the sample level.
            """,
            plot=table.plot(
                table_data,
                pconfig={
                    "id": "dragen-gc-metrics-summary-table",
                    "namespace": DragenGcMetrics.NAMESPACE,
                },
            ),
        )

        return data_by_sample.keys()

    @staticmethod
    def __get_normalized_gc_data(data_by_sample) -> dict:
        """Returns headers, data"""
        analysis = "GC BIAS DETAILS"
        hist_data = {}
        for sample_name, sample_data in data_by_sample.items():
            # {Normalized coverage at GC 0: 0.8324, Normalized coverage at GC 1,0.9456}
            hist_data[sample_name] = {
                int(key.split()[-1]): item
                for key, item in sample_data[analysis].items()
                if key.startswith("Normalized coverage at GC")
            }
        return hist_data

    @staticmethod
    def __get_summary_gc_data(data_by_sample) -> dict:
        analysis = "GC METRICS SUMMARY"
        summary_data = {}
        for sample_name, sample_data in data_by_sample.items():
            summary_data[sample_name] = {metric: stat for metric, stat in sample_data[analysis].items()}

        return summary_data


def parse_gc_metrics_file(f):
    """
    sample.gc_metrics.csv

    GC BIAS DETAILS,,Windows at GC 0,20,0.001
    GC BIAS DETAILS,,Windows at GC 1,4,0.000
    ...
    GC BIAS DETAILS,,Normalized coverage at GC 0,0.8324
    GC BIAS DETAILS,,Normalized coverage at GC 1,0.8540
    GC BIAS DETAILS,,Normalized coverage at GC 2,0.8512
    ...
    GC METRICS SUMMARY,,Window size,100
    GC METRICS SUMMARY,,Number of valid windows,3118404
    GC METRICS SUMMARY,,Number of discarded windows,366533633
    GC METRICS SUMMARY,,Average reference GC,40.84
    GC METRICS SUMMARY,,Mean global coverage,30.15
    GC METRICS SUMMARY,,Normalized coverage at GCs 0-19,1.06
    GC METRICS SUMMARY,,Normalized coverage at GCs 20-39,1.05
    GC METRICS SUMMARY,,Normalized coverage at GCs 40-59,0.97
    GC METRICS SUMMARY,,Normalized coverage at GCs 60-79,0.84
    GC METRICS SUMMARY,,Normalized coverage at GCs 80-100,0.47
    GC METRICS SUMMARY,,AT Dropout,0.58
    GC METRICS SUMMARY,,GC Dropout,2.01
    """

    data = defaultdict(dict)
    for line in f["f"].splitlines():
        tokens = line.split(",")
        # Percentage is currently unused
        if len(tokens) == 4:
            analysis, _, metric, stat = tokens
            percentage = None
        elif len(tokens) == 5:
            analysis, _, metric, stat, percentage = tokens
        else:
            raise ValueError(f"Unexpected number of tokens in line {line}")

        try:
            stat = float(stat)
        except ValueError:
            pass
        data[analysis][metric] = stat

    return data
