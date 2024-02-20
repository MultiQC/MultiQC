import logging

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.dragen.utils import Metric, make_headers
from multiqc.plots import table

log = logging.getLogger(__name__)

METRIC_NAMES = [
    "Unique cell-barcodes",
    "UMI threshold for passing cells",
    "Passing cells",
    "Median reads per cell",
    "Median UMIs per cell",
    "Median genes per cell",
    "Total genes detected",
]
METRICS = [
    Metric(
        id=m,
        title=" ".join([w.title() if w.islower() else w for w in m.split()]),
        in_genstats="#",
        in_own_tabl="#",
        precision=0,
        descr=m,
    )
    for m in METRIC_NAMES
]


class DragenScRnaMetrics(BaseMultiqcModule):
    def add_sc_rna_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/sc_rna_metrics"):
            data = parse_scrna_metrics_file(f)
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="sc_rna_metrics")
            data_by_sample[s_name] = data

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        # TODO: No test data for this file. No idea if we need this or not.
        # self.write_data_file(data_by_sample, "dragen_sc_rna_metrics")

        gen_stats_headers, table_headers = make_headers(METRIC_NAMES, METRICS)

        self.general_stats_addcols(data_by_sample, gen_stats_headers)
        self.add_section(
            name="Single-Cell RNA Metrics",
            anchor="sc-rna-metrics",
            description="""
            Summary metrics for single-cell RNA.
            """,
            plot=table.plot(
                {
                    sample_name: {
                        metric_name: int(stat) for metric_name, stat in metric.items() if metric_name in METRIC_NAMES
                    }
                    for sample_name, metric in data_by_sample.items()
                },
                table_headers,
                pconfig={
                    "namespace": "Single-Cell RNA Metrics",
                    "id": "dragen-sc-rna-metrics-table",
                },
            ),
        )

        return data_by_sample.keys()


def parse_scrna_metrics_file(f):
    """
    sample.scRNA.metrics.csv

    RUN TIME,,Time loading reference,00:01:31.289,91.29
    RUN TIME,,Time aligning reads,00:00:25.190,25.19
    RUN TIME,,Time duplicate marking,00:00:01.817,1.82
    RUN TIME,,Time sorting and marking duplicates,00:00:07.368,7.37
    RUN TIME,,Time DRAGStr calibration,00:00:07.069,7.07
    """
    data = {}
    for line in f["f"].splitlines():
        tokens = line.split(",")
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
        data[metric] = stat

    return data
