import logging

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
log = logging.getLogger(__name__)


NAMESPACE = "Ploidy estimation"


class DragenPloidyEstimationMetrics(BaseMultiqcModule):
    def add_ploidy_estimation_metrics(self):
        data_by_sample = dict()

        for f in self.find_log_files("dragen/ploidy_estimation_metrics"):
            data = parse_ploidy_estimation_metrics_file(f)
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            self.add_data_source(f, section="ploidy_estimation_metrics")
            data_by_sample[s_name] = data

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_sample, "dragen_ploidy")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        headers = {
            "Ploidy estimation": {
                "title": "Sex",
                "description": "Sex chromosome ploidy estimation (XX, XY, X0, 00, etc.)",
                "scale": "Set3",
            }
        }
        self.general_stats_addcols(data_by_sample, headers, namespace=NAMESPACE)

        table_headers = dict(
            **headers,
            **{
                "Autosomal median coverage": {
                    "title": "Autosomal med. cov.",
                    "description": "Median coverage of autosomal chromosomes",
                    "scale": "Blues",
                    "format": "{:.2f}x",
                },
                "X median coverage": {
                    "title": "X med. cov.",
                    "description": "Median coverage of X chromosome",
                    "scale": "Blues",
                    "format": "{:.2f}x",
                },
                "Y median coverage": {
                    "title": "Y med. cov.",
                    "description": "Median coverage of Y chromosome",
                    "scale": "Blues",
                    "format": "{:.2f}x",
                },
                "X median / Autosomal median": {
                    "title": "X med. / autosomal med. cov.",
                    "description": "Ratio of median coverage of X chromosome to median coverage of autosomal chromosomes",
                    "scale": "Blues",
                    "format": "{:.2f}",
                },
                "Y median / Autosomal median": {
                    "title": "Y med. / autosomal med. cov.",
                    "description": "Ratio of median coverage of Y chromosome to median coverage of autosomal chromosomes",
                    "scale": "Blues",
                    "format": "{:.2f}",
                },
            },
        )

        self.add_section(
            name="Ploidy estimation metrics",
            anchor="dragen_ploidy",
            plot=table.plot(
                data_by_sample,
                headers=table_headers,
                pconfig={
                    "id": "dragen_ploidy_table",
                    "title": "Ploidy estimation metrics",
                    "namespace": NAMESPACE,
                },
            ),
        )

        return data_by_sample.keys()


def parse_ploidy_estimation_metrics_file(f):
    """
    T_SRR7890936_50pc.ploidy_estimation_metrics.csv

    PLOIDY ESTIMATION,,Autosomal median coverage,55.63
    PLOIDY ESTIMATION,,X median coverage,27.44
    PLOIDY ESTIMATION,,Y median coverage,0.00
    PLOIDY ESTIMATION,,X median / Autosomal median,0.49
    PLOIDY ESTIMATION,,Y median / Autosomal median,0.00
    PLOIDY ESTIMATION,,Ploidy estimation,X0
    """

    data = {}
    for line in f["f"].splitlines():
        _, _, metric, stat = line.split(",")
        if stat.strip() == "":
            continue
        try:
            stat = float(stat)
        except ValueError:
            pass
        data[metric] = stat

    return data
