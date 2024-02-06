""" MultiQC submodule to parse output from deepTools plotEnrichment """

import logging

from multiqc.plots import linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class PlotEnrichmentMixin:
    def parse_plot_enrichment(self):
        """Find plotEnrichment output."""
        self.deeptools_plotEnrichment = dict()
        for f in self.find_log_files("deeptools/plotEnrichment"):
            parsed_data = self.parsePlotEnrichment(f)
            for k, v in parsed_data.items():
                if k in self.deeptools_plotEnrichment:
                    log.warning(f"Replacing duplicate sample {k}.")
                self.deeptools_plotEnrichment[k] = v

            if len(parsed_data) > 0:
                self.add_data_source(f, section="plotEnrichment")

        self.deeptools_plotEnrichment = self.ignore_samples(self.deeptools_plotEnrichment)

        if len(self.deeptools_plotEnrichment) == 0:
            return 0

        # Write data to file
        self.write_data_file(self.deeptools_plotEnrichment, "deeptools_plot_enrich")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        dCounts = dict()
        dPercents = dict()
        for sample, v in self.deeptools_plotEnrichment.items():
            dCounts[sample] = dict()
            dPercents[sample] = dict()
            for category, v2 in v.items():
                dCounts[sample][category] = v2["count"]
                dPercents[sample][category] = v2["percent"]
        config = {
            "data_labels": [
                {"name": "Counts in features", "ylab": "Counts in feature"},
                {"name": "Percents in features", "ylab": "Percent of reads in feature"},
            ],
            "id": "deeptools_enrichment_plot",
            "title": "deepTools: Signal enrichment per feature",
            "ylab": "Counts in feature",
            "categories": True,
            "ymin": 0.0,
        }
        self.add_section(
            name="Feature enrichment",
            description="Signal enrichment per feature according to plotEnrichment",
            anchor="deeptools_enrichment",
            plot=linegraph.plot([dCounts, dPercents], pconfig=config),
        )

        return len(self.deeptools_plotEnrichment)

    def parsePlotEnrichment(self, f):
        d = {}
        firstLine = True
        for line in f["f"].splitlines():
            if firstLine:
                firstLine = False
                continue
            cols = line.strip().split("\t")

            if len(cols) != 5:
                log.warning(
                    "{} was initially flagged as the output from plotEnrichment, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()

            s_name = self.clean_s_name(cols[0], f)
            if s_name not in d:
                d[s_name] = dict()
            cols[1] = str(cols[1])
            if cols[1] in d[s_name]:
                log.warning(f"Replacing duplicate sample:featureType {s_name}:{cols[1]}.")
            d[s_name][cols[1]] = dict()

            try:
                d[s_name][cols[1]]["percent"] = float(cols[2])
                d[s_name][cols[1]]["count"] = self._int(cols[3])
            except:  # noqa: E722
                log.warning(
                    "{} was initially flagged as the output from plotEnrichment, but that seems to not be the case. Skipping...".format(
                        f["fn"]
                    )
                )
                return dict()
        return d
