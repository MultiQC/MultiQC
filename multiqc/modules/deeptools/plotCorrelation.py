"""MultiQC submodule to parse output from deepTools plotCorrelation"""

import logging

from multiqc.plots import heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class plotCorrelationMixin:
    def parse_plotCorrelation(self):
        """Find plotCorrelation output"""
        self.deeptools_plotCorrelationData = dict()
        for f in self.find_log_files("deeptools/plotCorrelationData", filehandles=False):
            parsed_data = self.parsePlotCorrelationData(f)
            for sample, val_by_sample in parsed_data.items():
                if sample in self.deeptools_plotCorrelationData:
                    log.warning(f"Replacing duplicate sample {sample}.")
                self.deeptools_plotCorrelationData[sample] = val_by_sample
            if len(parsed_data) > 0:
                self.add_data_source(f, section="plotCorrelation")

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        self.deeptools_plotCorrelationData = self.ignore_samples(self.deeptools_plotCorrelationData)
        for s_name, val_by_sample in self.deeptools_plotCorrelationData.items():
            self.deeptools_plotCorrelationData[s_name] = self.ignore_samples(val_by_sample)

        if len(self.deeptools_plotCorrelationData) > 0:
            # Data is in wrong format for writing to file
            # self.write_data_file(self.deeptools_plotCorrelationData, "deeptools_plot_corr")

            config = {
                "id": "deeptools_correlation_plot",
                "title": "deeptools: Correlation Plot",
            }
            samples = sorted(self.deeptools_plotCorrelationData.keys())
            rows = []
            for s1 in samples:
                row = []
                for s2 in samples:
                    try:
                        row.append(self.deeptools_plotCorrelationData[s1][s2])
                    except KeyError:
                        row.append(None)
                rows.append(row)
            if len(rows) == 0:
                log.debug("No valid data for correlation plot")
                return None

            self.add_section(
                name="Correlation heatmap",
                anchor="deeptools_correlation",
                description="Pairwise correlations of samples based on distribution of sequence reads",
                plot=heatmap.plot(rows, samples, samples, config),
            )

        return len(self.deeptools_plotCorrelationData)

    def parsePlotCorrelationData(self, f):
        d = dict()
        x_samples = None
        for line in f["f"].splitlines():
            cols = line.split("\t")
            if cols[0] == "#plotCorrelation --outFileCorMatrix":
                continue
            elif cols[0] == "":
                x_samples = [self.clean_s_name(s.strip("'")) for s in cols[1 : len(cols)]]
            else:
                y_sample = str(cols[0]).strip("'")
                y_sample = self.clean_s_name(y_sample, f)
                d[y_sample] = dict()
                for x_sample, col in zip(x_samples, cols[1 : len(cols)]):
                    d[y_sample][x_sample] = float(col)
        return d
