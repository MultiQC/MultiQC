import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class DragenTimeMetrics(BaseMultiqcModule):
    def add_time_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/time_metrics"):
            data = parse_time_metrics_file(f)
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="time_metrics")
            data_by_sample[s_name] = data

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)

        if not data_by_sample:
            return set()

        # Save data file
        self.write_data_file(data_by_sample, "dragen_time_metrics")

        self.add_section(
            name="Time Metrics",
            anchor="dragen-time-metrics",
            description="Time metrics for DRAGEN run.  Total run time is less than the sum of individual "
            "steps because of parallelization.",
            plot=bargraph.plot(
                [
                    {
                        sample: {
                            step_name: step_time / 60
                            for step_name, step_time in data.items()
                            if step_name == "Total runtime"
                        }
                        for sample, data in data_by_sample.items()
                    },
                    {
                        sample: {
                            step_name[5:]: step_time / 60
                            for step_name, step_time in data.items()
                            if step_name != "Total runtime"
                        }
                        for sample, data in data_by_sample.items()
                    },
                ],
                pconfig={
                    "id": "time_metrics_plot",
                    "title": "Dragen: Time Metrics",
                    "ylab": "Time (minutes)",
                    "cpswitch_counts_label": "Time (minutes)",
                    "data_labels": [
                        {
                            "name": "Total Runtime",
                            "ylab": "Time (minutes)",
                            "cpswitch_counts_label": "Time (minutes)",
                        },
                        {
                            "name": "Steps Breakdown",
                            "ylab": "Time (minutes)",
                            "cpswitch_counts_label": "Time (minutes)",
                        },
                    ],
                },
            ),
        )

        return data_by_sample.keys()


def parse_time_metrics_file(f):
    """
    sample.time_metrics.csv

    RUN TIME,,Time loading reference,00:01:31.289,91.29
    RUN TIME,,Time aligning reads,00:00:25.190,25.19
    RUN TIME,,Time duplicate marking,00:00:01.817,1.82
    RUN TIME,,Time sorting and marking duplicates,00:00:07.368,7.37
    RUN TIME,,Time DRAGStr calibration,00:00:07.069,7.07
    """
    data = defaultdict(dict)
    for line in f["f"].splitlines():
        tokens = line.split(",")
        analysis, _, metric, timestr, seconds = tokens

        try:
            seconds = float(seconds)
        except ValueError:
            pass
        data[metric] = seconds

    return data
