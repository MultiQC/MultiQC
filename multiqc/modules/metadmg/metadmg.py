""" MultiQC module to plot output from metadmg """


import logging
import pandas as pd

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    To allow reading gzip archives, run with `ignore_images: false`
    in the config, e.g.:
    ```
    multiqc . --cl-config 'ignore_images: false'
    ```
    """

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="metaDMG",
            anchor="metadmg",
            href="https://github.com/metaDMG-dev/metaDMG-cpp",
            info="Taxonomic classification of reads and estimation of damage rates in ancient DNA data",
            doi="10.1101/2022.12.06.519264",
        )

        # Read config
        self.rank = getattr(config, "metadmg", {}).get("rank", "genus")
        self.n_taxa = getattr(config, "metadmg", {}).get("n_taxa", 20)
        self.sort_by = getattr(config, "metadmg", {}).get("sort_by", "nreads")
        self.index = getattr(config, "metadmg", {}).get("index", "name")

        # Read files
        self.metadmg_data = self.metadmg_read_data()

        # Filter to strip out ignored sample names
        self.metadmg_data = self.ignore_samples(self.metadmg_data)

        if len(self.metadmg_data) == 0:
            raise ModuleNoSamplesFound

        # Write parsed report data to a file
        self.write_data_file(self.metadmg_data, "multiqc_metadmg")

        # Add version
        self.add_software_version(None)

        # Plots section
        self.barplot_section()
        #self.dfitplot_section()

    def metadmg_read_data(self):
        metadmg_df = dict()

        for f in self.find_log_files("metadmg/stat", filecontents=False):
            f["s_name"] = f["s_name"][:-5]
            # Read DF
            metadmg_df[f["s_name"]] = pd.read_table(
                f["fn"],
                sep="\t",
                index_col="taxid",
                comment="#",
            )
            self.add_data_source(f, f["s_name"])

        for f in self.find_log_files("metadmg/dfit", filecontents=False):
            f["s_name"] = f["s_name"][:-5]
            # Read DF
            metadmg_df[f["s_name"]] = metadmg_df[f["s_name"]].join(
                pd.read_table(
                    f["fn"],
                    sep="\t",
                    index_col="taxid",
                    comment="#",
                )
            )
            self.add_data_source(f, f["s_name"])

        metadmg_data = dict()
        for s_name, df in metadmg_df.items():
            # Filter DF on rank
            df = df[df["rank"] == self.rank].drop(columns="rank")
            # Get first N taxa
            df = df.head(n=self.n_taxa).set_index("name")

            parsed_data = df.to_dict()
            if len(parsed_data) > 0:
                if s_name in metadmg_data:
                    log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
                metadmg_data[s_name] = parsed_data

        return metadmg_data

    def barplot_section(self):
        # Convert data
        stat_labels = {
            "nreads": "Number of reads",
            "nalign": "Number of alignments",
            "mean_rlen": "Mean read length",
            "var_rlen": "Read length variance",
            "mean_gc": "Mean GC content",
            "var_gc": "GC content variance",
            "A": "Damage estimate",
            "A_b": "Damage estimate (bootstrap)",
        }
        
        data_plot = list()
        data_labels = list()
        for stat_type in stat_labels.keys():
            data_plot.append({s_name: data[stat_type] for s_name, data in self.metadmg_data.items() if stat_type in data})
            data_labels.append({"name": stat_labels[stat_type]})

        # Config for the plot
        pconfig = {
            "id": "metadmg_rank_plot",
            "hide_zero_cats": False,
            "title": "metaDMG: read statistics",
            "ylab": None,
            "use_legend": False,
            "tt_decimals": 2,
            "tt_percentages": False,
            "data_labels": data_labels,
        }

        self.add_section(
            name="Read statistics by taxonomic rank",
            anchor="metadmg-rank",
            description=f"Read abundance statistics for top {self.n_taxa} {self.rank}",
            plot=bargraph.plot(data_plot, pconfig=pconfig),
        )

    def dfitplot_section(self):
	# Convert data
        data_plot = list()
        data_labels = list()
        for s_name, data in sorted(self.metadmg_data.items()):
            log.debug(data)
            data_plot.append(data)
            data_labels.append({"name": s_name})

        # Config for the plot
        pconfig = {
            "id": "metadmg_dfit_plot",
            "title": "metaDMG: damage estimation",
            "ylab": None,
            "data_labels": data_labels,
        }

        self.add_section(
            name="Damage estimates by taxonomic rank",
            anchor="metadmg-dfit",
            description=f"Damage estimates for top {self.n_taxa} {self.rank}",
            plot=linegraph.plot(data_plot, pconfig=pconfig),
        )
