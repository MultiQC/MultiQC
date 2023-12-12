""" MultiQC module to parse output from nonpareil """


import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# Initialise the logger
from multiqc.plots import linegraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Nonpareil analysis is split into two parts: the first (written in C++)
    performs the subsampling, while the second (written in R) performs the
    statistical analyses. As such, this model requires the user to post-process
    the R object from the second part."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="nonpareil",
            anchor="nonpareil",
            href="https://github.com/lmrodriguezr/nonpareil",
            info="Estimate metagenomic coverage and sequence diversity ",
            doi="10.1093/bioinformatics/btt584",
        )
        # Config options
        self.plot_observed = getattr(config, "nonpareil", {}).get("plot_observed", True)
        self.plot_model = getattr(config, "nonpareil", {}).get("plot_model", True)
        self.disp_type = getattr(config, "nonpareil", {}).get("plot_dispersion", False)

        # Read JSON file
        self.data_by_sample = dict()
        for f in self.find_log_files("nonpareil", filehandles=True):
            json_parsed = self.parse_nonpareil_json(f)
            common_samples = set(json_parsed).intersection(self.data_by_sample)
            if common_samples:
                log.debug(f"Duplicate sample names found ({common_samples})! Overwriting...")
            self.data_by_sample = self.data_by_sample | json_parsed
        # Register samples
        for s_name in self.data_by_sample.keys():
            self.add_data_source(f, s_name)

        # Filter to strip out ignored sample names
        self.data_by_sample = self.ignore_samples(self.data_by_sample)
        if len(self.data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.data_by_sample)} reports")
        self.write_data_file(self.data_by_sample, "nonpareil")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        version = set([data["version"] for sample, data in self.data_by_sample.items()])
        if len(version) != 1:
            raise ValueError("several versions found.")
        self.add_software_version(str(version.pop()))

        # Add general stats
        self.nonpareil_general_stats_table()

        # Add plots
        self.nonpareil_redundancy_plot()

    def parse_nonpareil_json(self, f):
        """Parse the JSON output from nonpareil"""
        import json

        try:
            json_raw = json.load(f["f"])
        except json.JSONDecodeError as e:
            log.warning(f"Could not parse nonpareil JSON: '{f['fn']}': {e}, skipping file")
            return None, {}

        # x.obs Rarefied sequencing effort.
        # x.adj Adjusted rarefied sequencing effort.
        # y.red Rarefied redundancy (observed).
        # y.cov Rarefied coverage (corrected).
        # y.sd Standard deviation of rarefied coverage.
        # y.p25 Percentile 25 (1st quartile) of rarefied coverage.
        # y.p50 Percentile 50 (median) of rarefied coverage.
        # y.p75 Percentile 75 (3rd quartile) of rarefied coverage

        import numpy as np

        for s_name, data in json_raw.items():
            # Conbert base pairs to gigabase pairs
            data["LR"] = data["LR"] / 1000 / 1000
            data["LRstar"] = data["LRstar"] / 1000 / 1000
            data["x.adj"] = [(1e-6 if x == 0 else x) / 1000 / 1000 for x in data["x.adj"]]
            data["x.model"] = [x / 1000 / 1000 for x in data["x.model"]]
            data["observed"] = {x: y for x, y in zip(data["x.adj"], data["y.cov"])}
            data["model"] = {x: y for x, y in zip(data["x.model"], data["y.model"])}
            # Calculate dispersion
            # from https://github.com/lmrodriguezr/nonpareil/blob/162f1697ab1a21128e1857dd87fa93011e30c1ba/utils/Nonpareil/R/Nonpareil.R#L306-L318
            disp_add = False
            if self.disp_type == "sd":
                disp_add = data["y.sd"]
            elif self.disp_type == "ci95":
                disp_add = list(np.array(data["y.sd"]) * 1.9)
            elif self.disp_type == "ci90":
                disp_add = list(np.array(data["y.sd"]) * 1.64)
            elif self.disp_type == "ci50":
                disp_add = list(np.array(data["y.sd"]) * 0.67)
            elif self.disp_type == "iq":
                data[f"{self.disp_type}_upper"] = {x: y for x, y in zip(data["x.adj"], data["y.p75"])}
                data[f"{self.disp_type}_lower"] = {x: y for x, y in zip(data["x.adj"], data["y.p25"])}

            if disp_add:
                data[f"{self.disp_type}_upper"] = {
                    x: y for x, y in zip(data["x.adj"], list(np.array(data["y.cov"]) + np.array(disp_add)))
                }
                data[f"{self.disp_type}_lower"] = {
                    x: y for x, y in zip(data["x.adj"], list(np.array(data["y.cov"]) - np.array(disp_add)))
                }

        return json_raw

    def nonpareil_general_stats_table(self):
        """Take the parsed stats from the nonpareil report and add it to the
        General Statistics table at the top of the report"""

        headers = {
            "L": {
                "title": "Read Length",
                "min": 0,
                "suffix": " bps",
                "scale": "GnBu",
                "format": "{:,.3f}",
                "hidden": True,
            },
            "AL": {
                "title": "Adjusted read length",
                "description": "Adjusted read length (if using kmers)",
                "min": 0,
                "scale": "GnBu",
                "suffix": " bps",
                "format": "{:,.3f}",
                "hidden": True,
            },
            "R": {
                "title": "Number of reads",
                "min": 0,
                "scale": "GnBu",
                "format": "{:,.0f}",
                "hidden": True,
            },
            "LR": {
                "title": "Sequencing effort",
                "description": "Effective sequencing effort used",
                "min": 0,
                "suffix": " Mbps",
                "scale": "Blues",
                "format": "{:,.2f}",
                "hidden": True,
            },
            "overlap": {
                "title": "Read overlap",
                "description": "Minimum read overlap",
                "min": 0,
                "suffix": " bps",
                "scale": "Blues",
                "hidden": True,
                "format": "{:,.0f}",
            },
            "log.sample": {
                "title": "Log sampling",
                "description": "Multiplier of the log-sampling (or zero if linear)",
                "max": 1,
                "min": 0,
                "scale": "BuGn",
                "hidden": True,
                "format": "{:,.2f}",
            },
            "kappa": {
                "title": "Redundancy",
                "description": "Dataset redundancy",
                "max": 1,
                "min": 0,
                "scale": "RdYlGn-rev",
                "format": "{:,.2f}",
            },
            "C": {
                "title": "Coverage",
                "description": "Dataset coverage",
                "max": 1,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.2f}",
            },
            "consistent": {
                "title": "Consistent Data?",
                "description": "Is the data sufficient for accurate estimation?",
                "modify": lambda x: "Yes" if x == 1 else "No",
                "scale": False,
                "hidden": True,
            },
            "star": {
                "title": "Ideal coverage",
                "description": "Coverage considered ’nearly complete’.",
                "min": 0,
                "suffix": "%",
                "scale": False,
                "hidden": True,
            },
            "LRstar": {
                "title": "Sequencing effort for ideal coverage",
                "description": "Projected sequencing effort for nearly complete coverage",
                "min": 0,
                "suffix": " Mbps",
                "scale": "RdYlGn-rev",
                "format": "{:,.2f}",
            },
            "modelR": {
                "title": "Model correlation",
                "description": "Pearson R for the estimated model",
                "max": 1,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.3f}",
                "hidden": True,
            },
            "diversity": {
                "title": "Diversity",
                "description": "Dataset Nd index of sequence diversity",
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.2f}",
            },
        }

        self.general_stats_addcols(self.data_by_sample, headers)

    def nonpareil_redundancy_plot(self):
        """Make the redundancy plot for nonpareil"""

        data_plot = list()
        data_labels = list()
        for s_name, data in sorted(self.data_by_sample.items()):
            data_plot.append(dict())
            if self.plot_observed:
                data_plot[-1]["observed"] = data["observed"]
            if self.plot_model:
                data_plot[-1]["model"] = data["model"]
            if self.disp_type:
                data_plot[-1][f"{self.disp_type}_upper"] = data[f"{self.disp_type}_upper"]
                data_plot[-1][f"{self.disp_type}_lower"] = data[f"{self.disp_type}_lower"]
            data_labels.append({"name": s_name})

        log.debug(data_labels)

        pconfig = {
            "id": "nonpareil-redundancy-plot",
            "title": "Nonpareil: Redundancy levels",
            "xlab": "Sequencing effort (Mbps)",
            "ylab": "Estimated Average Coverage",
            "xmin": 1e-3,
            "ymin": 0,
            "ymax": 1,
            "yCeiling": 1,
            "xDecimals": True,
            "yDecimals": True,
            "xLog": True,
            "tt_label": "{point.x:.2f} Mbps: {point.y:.2f}",
            "data_labels": data_labels,
        }
        log.debug(pconfig)

        self.add_section(
            name="Redundancy levels",
            anchor="nonpareil-redundancy",
            description="Redundancy levels across samples.",
            plot=linegraph.plot(data_plot, pconfig),
        )
