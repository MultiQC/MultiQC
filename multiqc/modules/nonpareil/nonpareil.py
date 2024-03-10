""" MultiQC module to parse output from nonpareil """


import logging
import numpy as np

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.utils import mqc_colour


# Initialise the logger
from multiqc.plots import table, linegraph

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
        self.plot_colours = getattr(config, "nonpareil", {}).get("plot_colours", "Paired")

        # Read JSON file
        self.data_by_sample = dict()
        for f in self.find_log_files("nonpareil", filehandles=True):
            json_parsed = self.parse_nonpareil_json(f)
            common_samples = set(json_parsed).intersection(self.data_by_sample)
            if common_samples:
                log.debug(f"Duplicate sample names found ({common_samples})! Overwriting...")
            self.data_by_sample.update(json_parsed)
            for s_name in json_parsed:
                self.add_data_source(f, s_name)

        # Filter to strip out ignored sample names
        self.data_by_sample = self.ignore_samples(self.data_by_sample)
        if len(self.data_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(self.data_by_sample)} reports")
        self.write_data_file(self.data_by_sample, "nonpareil")

        # Add versions
        for s_name, data in self.data_by_sample.items():
            self.add_software_version(str(data["nonpareil_version"]), s_name)

        # Add general stats
        self.nonpareil_general_stats_table()
        # Add table
        self.nonpareil_stats_table()
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

        # Rename samples, if label available
        for s_name in list(json_raw.keys()):
            s_label = json_raw[s_name].get("label")
            if s_label:
                json_raw[s_label] = json_raw.pop(s_name)

        for s_name, data in json_raw.items():
            # Convert base pairs to megabase pairs
            data["x.adj"] = [(1e-6 if x == 0 else x) / 1000 / 1000 for x in data["x.adj"]]
            if data["has.model"]:
                data["x.model"] = [x / 1000 / 1000 for x in data["x.model"]]
            # Convert fraction to percentage
            data["y.cov"] = [x * 100 for x in data["y.cov"]]
            if data["has.model"]:
                data["y.model"] = [x * 100 for x in data["y.model"]]
            data["y.sd"] = [x * 100 for x in data["y.sd"]]
            data["y.p75"] = [x * 100 for x in data["y.p75"]]
            data["y.p25"] = [x * 100 for x in data["y.p25"]]
            # Prepare plot data
            data["observed"] = {x: y for x, y in zip(data["x.adj"], data["y.cov"])}
            if data["has.model"]:
                data["model"] = {x: y for x, y in zip(data["x.model"], data["y.model"])}
            else:
                modelR = data["modelR"]
                assert isinstance(modelR, list) and len(modelR) == 0, "there is no model, but modelR is not empty"
                data["modelR"] = np.nan
                LRstar = data["LRstar"]
                assert isinstance(LRstar, list) and len(LRstar) == 0, "there is no model, but LRstar is not empty"
                data["LRstar"] = np.nan
            # Add prefix to labels
            for key in list(data.keys()):
                data["nonpareil_" + key] = data.pop(key)

        return json_raw

    def nonpareil_general_stats_table(self):
        """Take the parsed stats from the nonpareil report and add it to the
        General Statistics table at the top of the report"""

        headers = {
            "nonpareil_R": {
                "title": f"{config.read_count_prefix} Reads",
                "description": f"Total raw sequences ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "min": 0,
                "scale": "RdYlGn",
                "shared_key": "read_count",
                "hidden": True,
            },
            "nonpareil_LR": {
                "title": f"{config.base_count_prefix} Seq. effort",
                "description": f"Total base pairs sequenced ({config.base_count_desc})",
                "modify": lambda x: x * config.base_count_multiplier,
                "min": 0,
                "scale": "Blues",
                "shared_key": "base_count",
                "hidden": True,
            },
            "nonpareil_kappa": {
                "title": "Redundancy",
                "description": "Dataset redundancy",
                "modify": lambda x: x * 100,
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
            "nonpareil_C": {
                "title": "Coverage",
                "description": "Dataset coverage",
                "modify": lambda x: x * 100,
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "nonpareil_diversity": {
                "title": "Diversity",
                "description": "Dataset Nd index of sequence diversity",
                "min": 0,
                "scale": "GnBu-rev",
            },
        }

        self.general_stats_addcols(self.data_by_sample, headers)

    def nonpareil_stats_table(self):
        """Take the parsed stats from the nonpareil report and add it to the
        General Statistics table at the top of the report"""

        headers = {
            "nonpareil_L": {
                "title": "Read Length",
                "min": 0,
                "suffix": " bp",
                "scale": "GnBu",
                "format": "{:,.3f}",
                "hidden": True,
            },
            "nonpareil_AL": {
                "title": "Adjusted read length",
                "description": "Adjusted read length (if using kmers - same as read length for alignment)",
                "min": 0,
                "scale": "GnBu",
                "suffix": " bp",
                "format": "{:,.3f}",
                "hidden": True,
            },
            "nonpareil_R": {
                "title": f"{config.read_count_prefix} Reads",
                "description": f"Total raw sequences ({config.read_count_desc})",
                "modify": lambda x: x * config.read_count_multiplier,
                "min": 0,
                "scale": "RdYlGn",
                "shared_key": "read_count",
            },
            "nonpareil_LR": {
                "title": f"{config.base_count_prefix} Seq. effort",
                "description": f"Total base pairs sequenced ({config.base_count_desc})",
                "modify": lambda x: x * config.base_count_multiplier,
                "min": 0,
                "scale": "Blues",
                "shared_key": "base_count",
            },
            "nonpareil_overlap": {
                "title": "Read overlap",
                "description": "Minimum read overlap",
                "min": 0,
                "suffix": " bp",
                "hidden": True,
                "format": "{:,.0f}",
            },
            "nonpareil_log.sample": {
                "title": "Log sampling",
                "description": "Multiplier of the log-sampling (or zero if linear)",
                "max": 1,
                "min": 0,
                "hidden": True,
            },
            "nonpareil_kappa": {
                "title": "Redundancy",
                "description": "Dataset redundancy",
                "modify": lambda x: x * 100,
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn-rev",
            },
            "nonpareil_C": {
                "title": "Coverage",
                "description": "Dataset coverage",
                "modify": lambda x: x * 100,
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "nonpareil_consistent": {
                "title": "Consistent Data?",
                "description": "Is the data sufficient for accurate estimation?",
                "modify": lambda x: "Yes" if x == 1 else "No",
                "scale": False,
            },
            "nonpareil_star": {
                "title": "Ideal coverage",
                "description": "Objective coverage in percentage; i.e., coverage value considered near-complete.",
                "max": 100,
                "min": 0,
                "suffix": "%",
                "scale": False,
                "hidden": True,
            },
            "nonpareil_has.model": {
                "title": "Model available?",
                "description": "Was model estimated?",
                "modify": lambda x: "Yes" if x else "No",
                "scale": False,
                "hidden": True,
            },
            "nonpareil_LRstar": {
                "title": f"{config.base_count_prefix} Ideal seq. effort",
                "description": f"Projected sequencing effort for nearly complete coverage ({config.base_count_desc})",
                "modify": lambda x: x * config.base_count_multiplier,
                "min": 0,
                "scale": "RdYlGn-rev",
                "shared_key": "base_count",
            },
            "nonpareil_modelR": {
                "title": "Model correlation",
                "description": "Pearson R for the estimated model",
                "max": 1,
                "min": 0,
                "scale": "RdYlGn",
                "format": "{:,.3f}",
                "hidden": True,
            },
            "nonpareil_diversity": {
                "title": "Diversity",
                "description": "Dataset Nd index of sequence diversity",
                "min": 0,
                "scale": "GnBu-rev",
            },
        }

        pconfig = {
            "id": "nonpareil-table",
            "title": "Nonpareil: statistics",
        }

        self.add_section(
            name="Statistics",
            anchor="nonpareil-table-section",
            description="""
            Nonpareil uses the redundancy of the reads in metagenomic datasets to
            estimate the average coverage and predict the amount of sequences that
            will be required to achieve "nearly complete coverage".
            """,
            helptext="""
            From the Nonpareil [FAQ](http://enve-omics.ce.gatech.edu/nonpareil/faq):

            **What is the Nonpareil diversity index?**

            We haven't yet published a formal description of the "diversity" reported by Nonpareil, which we internally refer to as Nonpareil diversity. The Nonpareil diversity is a measurement of how complex a community is in "sequence space". Graphically, the index is a measurement of "how much to the right the Nonpareil curve is".

            **What's the expected range of the Nonpareil diversity index?**

            We've found that the Nonpareil diversity index ranks communities very
            reliably in terms of diversity, and can even detect small seasonal
            variations in bacterial communities.
            For reference: an environment largely dominated by a single bacterial
            species like human posterior fornix or the Acid Mine Drainage gets
            values around 15 to 16, human stool samples around 17 to 19, freshwater
            and sandy soils around 20 to 22, and marine and other soils around
            21 to 25. This index has a logarithmic scale, to allow meaningful
            comparisons between extreme cases, like single species vs soil.
            """,
            plot=table.plot(self.data_by_sample, headers, pconfig),
        )

    def nonpareil_redundancy_plot(self):
        """Make the redundancy plot for nonpareil"""

        extra_series_config = {
            "dash": "dash",
            "line": {"width": 2},
            "showlegend": False,
        }

        data_colors_default = mqc_colour.mqc_colour_scale().get_colours(self.plot_colours)
        data_colors = {
            s_name: data.get("nonpareil_col", data_colors_default.pop(0))
            for s_name, data in self.data_by_sample.items()
        }

        data_labels = [
            {"name": "Combined"},
            {"name": "Observed"},
            {"name": "Model"},
        ]
        data_plot = list()
        extra_series = list()
        for idx, dataset in enumerate(data_labels):
            data_plot.append(dict())
            extra_series.append(list())
            for s_name, data in self.data_by_sample.items():
                if dataset["name"] == "Observed":
                    data_plot[idx][s_name] = data["nonpareil_observed"]
                elif dataset["name"] == "Model" and data["nonpareil_has.model"]:
                    data_plot[idx][s_name] = data["nonpareil_model"]
                elif dataset["name"] == "Combined":
                    data_plot[idx][s_name] = data["nonpareil_observed"]
                    if data["nonpareil_has.model"]:
                        extra_series[idx].append(dict(extra_series_config))
                        extra_series[idx][-1]["name"] = s_name
                        extra_series[idx][-1]["data"] = [[x, y] for x, y in data["nonpareil_model"].items()]
                        extra_series[idx][-1]["color"] = data_colors[s_name]

        pconfig = {
            "id": "nonpareil-redundancy-plot",
            "colors": data_colors,
            "title": "Nonpareil: Redundancy levels",
            "xlab": "Sequencing effort (Mbp)",
            "ylab": "Estimated Average Coverage",
            "xmin": 1e-3,
            "ymin": 0,
            "ymax": 100,
            "xDecimals": True,
            "yDecimals": True,
            "xLog": True,
            "tt_label": "{point.x:.2f} Mbp: {point.y:.2f}",
            "extra_series": extra_series,
            "data_labels": data_labels,
        }

        self.add_section(
            name="Redundancy levels",
            anchor="nonpareil-redundancy",
            description="Observed and modelled redundancy levels across samples.",
            helptext="""
            The estimation of the Redundancy is at the core of Nonpareil,
            but it's when those values are transformed into average coverage
            that they become comparable across samples, and become useful for
            project design and sample evaluation.
            """,
            plot=linegraph.plot(data_plot, pconfig),
        )
