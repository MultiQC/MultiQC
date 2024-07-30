"""MultiQC submodule to parse output from deepTools plotFingerprint"""

import logging

import numpy as np

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class PlotFingerprintMixin:
    def parse_plotFingerprint(self):
        """Find plotFingerprint output. Both --outQualityMetrics and --outRawCounts"""
        self.deeptools_plotFingerprintOutQualityMetrics = dict()
        for f in self.find_log_files("deeptools/plotFingerprintOutQualityMetrics"):
            parsed_data = self.parsePlotFingerprintOutQualityMetrics(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotFingerprintOutQualityMetrics:
                    log.warning(f"Replacing duplicate sample {k}.")
                # Values are fractions - convert to percentages for consistency with other MultiQC output
                self.deeptools_plotFingerprintOutQualityMetrics[k] = {i: float(j) * 100.0 for i, j in v.items()}

            if len(parsed_data) > 0:
                self.add_data_source(f, section="plotFingerprint")

        self.deeptools_plotFingerprintOutRawCounts = dict()
        for f in self.find_log_files("deeptools/plotFingerprintOutRawCounts"):
            parsed_data = self.parsePlotFingerprintOutRawCounts(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotFingerprintOutRawCounts:
                    log.warning(f"Replacing duplicate sample {k}.")
                self.deeptools_plotFingerprintOutRawCounts[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section="plotFingerprint")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.deeptools_plotFingerprintOutRawCounts = self.ignore_samples(self.deeptools_plotFingerprintOutRawCounts)
        self.deeptools_plotFingerprintOutQualityMetrics = self.ignore_samples(
            self.deeptools_plotFingerprintOutQualityMetrics
        )

        if len(self.deeptools_plotFingerprintOutRawCounts) > 0:
            # Write data to file
            self.write_data_file(self.deeptools_plotFingerprintOutRawCounts, "deeptools_plot_fingerprint_counts")

            self.add_section(
                name="Fingerprint plot",
                anchor="deeptools_fingerprint",
                description="Signal fingerprint according to plotFingerprint",
                plot=linegraph.plot(
                    self.deeptools_plotFingerprintOutRawCounts,
                    {
                        "id": "deeptools_fingerprint_plot",
                        "title": "deepTools: Fingerprint plot",
                        "xmin": 0.0,
                        "xmax": 1.0,
                        "ymin": 0.0,
                        "ymax": 1.0,
                        "xlab": "rank",
                        "ylab": "Fraction w.r.t. bin with highest coverage",
                    },
                ),
            )

        if len(self.deeptools_plotFingerprintOutQualityMetrics) > 0:
            # Write data to file
            self.write_data_file(self.deeptools_plotFingerprintOutQualityMetrics, "deeptools_plot_fingerprint_metrics")

            self.add_section(
                name="Fingerprint quality metrics",
                anchor="plotFingerprint",
                description="Various quality metrics returned by plotFingerprint",
                plot=linegraph.plot(
                    self.deeptools_plotFingerprintOutQualityMetrics,
                    {
                        "id": "plotFingerprint_quality_metrics",
                        "title": "deepTools: Fingerprint quality metrics",
                        "ymin": 0,
                        "ymax": 100,
                        "ylab_format": "{value}%",
                        "ylab": "Percentage of fragments",
                        "categories": True,
                        "tt_label": "<strong>{point.x}</strong>: {point.y:.2f}%",
                    },
                ),
            )

        return len(self.deeptools_plotFingerprintOutQualityMetrics), len(self.deeptools_plotFingerprintOutRawCounts)

    def parsePlotFingerprintOutQualityMetrics(self, f):
        d = {}
        firstLine = True
        header = []
        for line in f["f"].splitlines():
            cols = line.strip().split("\t")

            if len(cols) < 7:
                log.warning(
                    "{} was initially flagged as the output from plotFingerprint --outQualityMetrics, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()

            if firstLine:
                header = [str(x) for x in cols[1:]]
                firstLine = False
                continue

            s_name = self.clean_s_name(cols[0], f)
            if s_name in d:
                log.warning(f"Replacing duplicate sample {s_name}.")
            d[s_name] = dict()

            try:
                for i, c in enumerate(cols[1:]):
                    if i >= len(header):
                        log.warning(
                            "{} was initially flagged as the output from plotFingerprint --outQualityMetrics, but that seems to not be the case. Skipping...".format(
                                f["fn"]
                            )
                        )
                        return dict()
                    if header[i] == "AUC" or header[i] == "Synthetic AUC":
                        continue
                    d[s_name][header[i]] = float(c)
            except:  # noqa: E722
                log.warning(
                    "{} was initially flagged as the output from plotFingerprint --outQualityMetrics, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()
        return d

    def parsePlotFingerprintOutRawCounts(self, f):
        d = dict()
        samples = []
        firstLine = True
        for line in f["f"].splitlines():
            cols = line.strip().split("\t")
            if cols[0] == "#plotFingerprint --outRawCounts":
                continue

            if firstLine:
                for c in cols:
                    c = str(c).strip("'")
                    s_name = self.clean_s_name(c, f)
                    d[s_name] = []
                    samples.append(s_name)
                firstLine = False
                continue

            for idx, c in enumerate(cols):
                d[samples[idx]].append(self._int(c))

        # Switch to numpy, get the normalized cumsum
        x = np.linspace(
            0, len(d[samples[0]]) - 1, endpoint=True, num=100, dtype=int
        )  # The indices into the vectors that we'll actually return for plotting
        xp = np.arange(len(d[samples[0]]) + 1) / float(len(d[samples[0]]) + 1)
        for k, v in d.items():
            v = np.array(v)
            v = np.sort(v)
            cs = np.cumsum(v)
            cs = cs / float(cs[-1])
            # Convert for plotting
            v2 = dict()
            v2[0.0] = 0.0
            for _ in x:
                v2[xp[_]] = cs[_]
            d[k] = v2
        return d
