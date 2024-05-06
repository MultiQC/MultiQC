from typing import Dict

import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

log = logging.getLogger(__name__)


class DragenTrimmerMetrics(BaseMultiqcModule):
    def add_trimmer_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/trimmer_metrics"):
            data = parse_trimmer_metrics_file(f)
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="trimmer_metrics")
            data_by_sample[s_name] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            return set()

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Save data file
        self.write_data_file(data_by_sample, "dragen_trimmer_metrics")

        table_data = DragenTrimmerMetrics.__get_table_data(data_by_sample)
        self.add_section(
            name="Trimmer Metrics",
            anchor="trimmer-metrics",
            description="""
            Metrics on trimmed reads.
            """,
            plot=table.plot(
                table_data,
                pconfig={
                    "id": "dragen-trimmer-metrics-table",
                    "namespace": "Trimmer Metrics",
                    "title": "DRAGEN: Trimmer Metrics",
                },
            ),
        )

        return data_by_sample.keys()

    @staticmethod
    def __get_table_data(data_by_sample: Dict) -> Dict:
        ANALYSIS = "TRIMMER STATISTICS"
        trimmer_data = {}
        for sample, data in data_by_sample.items():
            trimmer_data[sample] = {metric: stat for metric, stat in data[ANALYSIS].items()}
        return trimmer_data


def parse_trimmer_metrics_file(f):
    """
    sample.trimmer_metrics.csv

    TRIMMER STATISTICS,,Total input reads,18531840
    TRIMMER STATISTICS,,Total input bases,2779776000
    TRIMMER STATISTICS,,Total input bases R1,1389888000
    TRIMMER STATISTICS,,Total input bases R2,1389888000
    TRIMMER STATISTICS,,Average input read length,150
    TRIMMER STATISTICS,,Total trimmed reads,0,0.00
    TRIMMER STATISTICS,,Total trimmed bases,0,0.00
    """

    data = defaultdict(dict)
    for line in f["f"].splitlines():
        tokens = line.split(",")
        if len(tokens) == 4:
            analysis, _, metric, stat = tokens
            try:
                stat = int(stat)
            except ValueError:
                try:
                    stat = float(stat)
                except ValueError:
                    pass
            val = stat
        elif len(tokens) == 5:
            analysis, _, metric, stat, percentage = tokens
            try:
                percentage = float(percentage)
            except ValueError:
                pass
            val = percentage
        else:
            log.error(f"DRAGEN trimming metrics: unexpected number of tokens in line {line}")
            continue
        data[analysis][metric] = val

    return data
