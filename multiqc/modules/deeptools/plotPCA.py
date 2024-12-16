"""MultiQC submodule to parse output from deepTools plotPCA"""

import logging

from multiqc.plots import scatter

# Initialise the logger
log = logging.getLogger(__name__)


class plotPCAMixin:
    def parse_plotPCA(self):
        """Find plotPCA output"""
        self.deeptools_plotPCAData = dict()
        for f in self.find_log_files("deeptools/plotPCAData", filehandles=False):
            parsed_data = self.parsePlotPCAData(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotPCAData:
                    log.warning(f"Replacing duplicate sample {k}.")
                self.deeptools_plotPCAData[k] = v
            if len(parsed_data) > 0:
                self.add_data_source(f, section="plotPCA")

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, f["s_name"])

        self.deeptools_plotPCAData = self.ignore_samples(self.deeptools_plotPCAData)

        if len(self.deeptools_plotPCAData) > 0:
            # Write data to file
            self.write_data_file(self.deeptools_plotPCAData, "deeptools_plot_PCA")

            config = {
                "id": "deeptools_pca_plot",
                "title": "deeptools: PCA Plot",
                "xlab": "PC1",
                "ylab": "PC2",
                "tt_label": "PC1 {point.x:.2f}: PC2 {point.y:.2f}",
            }
            data = dict()
            for s_name in self.deeptools_plotPCAData:
                try:
                    data[s_name] = {
                        "x": self.deeptools_plotPCAData[s_name][1],
                        "y": self.deeptools_plotPCAData[s_name][2],
                    }
                except KeyError:
                    pass
            if len(data) == 0:
                log.debug("No valid data for PCA plot")
                return None

            self.add_section(
                name="PCA plot",
                anchor="deeptools_pca",
                description="PCA plot with the top two principal components calculated based on genome-wide distribution of sequence reads",
                plot=scatter.plot(data, config),
            )

        return len(self.deeptools_plotPCAData)

    def parsePlotPCAData(self, f):
        d = dict()
        samples = []
        for line in f["f"].splitlines():
            cols = line.strip().split("\t")
            if cols[0] == "#plotPCA --outFileNameData":
                continue
            elif cols[0] == "Component":
                for c in cols[1 : (len(cols) - 1)]:
                    c = str(c).strip("'")
                    s_name = self.clean_s_name(c, f)
                    d[s_name] = {}
                    samples.append(s_name)
            else:
                idx = 0
                compo = cols[0]
                for c in cols[1 : (len(cols) - 1)]:
                    d[samples[idx]][self._int(compo)] = float(c)
                    idx += 1
        return d
