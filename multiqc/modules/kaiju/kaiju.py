""" MultiQC module to parse output from Kaiju """

import logging

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="Kaiju",
            anchor="kaiju",
            href="http://kaiju.binf.ku.dk/",
            info="a fast and sensitive taxonomic classification for metagenomics.",
            doi="10.1038/ncomms11257",
        )

        # Set up data structures
        self.kaiju_data = {}
        self.taxonomic_ranks = ["species", "genus", "family", "order", "class", "phylum"]

        # Find and parse kaiju2table report
        for f in self.find_log_files("kaiju", filehandles=True):
            taxo_rank, parsed_data = self.parse_kaiju2table_report(f)
            if taxo_rank is not None and parsed_data is not None:
                # One report can have several samples, raise error if sample already viewed with this taxonomic rank
                s_names = parsed_data.keys()
                for s_name in s_names:
                    if taxo_rank in self.kaiju_data and s_name in self.kaiju_data[taxo_rank].keys():
                        log.debug(
                            "Duplicate sample found in logs at {} rank! Overwriting sample: {}".format(
                                taxo_rank, s_name
                            )
                        )

                    # Superfluous function call to confirm that it is used in this module
                    # Replace None with actual version if it is available
                    self.add_software_version(None, s_name)

                self.add_data_source(f)
                if taxo_rank in self.kaiju_data.keys():
                    self.kaiju_data[taxo_rank].update(parsed_data)
                else:
                    self.kaiju_data[taxo_rank] = parsed_data

        # Filters to strip out ignored sample names
        all_samples_name = []
        for taxo_rank in self.kaiju_data.keys():
            self.kaiju_data[taxo_rank] = self.ignore_samples(self.kaiju_data[taxo_rank])
            all_samples_name.extend(self.kaiju_data[taxo_rank].keys())
        self.all_s_name = set(all_samples_name)
        # Number of samples found
        try:
            num_samples = max([len(taxo_rank) for taxo_rank in self.kaiju_data.values()])
        except ValueError:
            # No log files so didn't get any taxo_ranks
            raise ModuleNoSamplesFound

        # no file found
        if num_samples == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {num_samples} reports")

        self.kaiju_total_pct = dict()
        self.kaiju_sample_total_readcounts = dict()
        self.kaiju_sample_unclassified = dict()

        self.sum_sample_counts()
        self.kaiju_stats_table()
        self.top_five_barplot()

    def parse_kaiju2table_report(self, f):
        """Search a kaiju with a set of regexes"""
        parsed_data = {}
        taxo_rank = None

        for line in f["f"]:
            if line.startswith("file\t"):
                continue
            (s_file, pct, reads, taxon_id, taxon_names) = line.rstrip().split("\t")
            s_name = self.clean_s_name(s_file, f)
            if s_name not in parsed_data.keys():
                parsed_data[s_name] = {"assigned": {}}

            if taxon_names.startswith("cannot be assigned") or taxon_names == "unclassified":
                if taxon_names.startswith("cannot be assigned"):
                    parsed_data[s_name]["cannot be assigned"] = {"count": int(reads), "percent": float(pct)}
                    taxo_rank = taxon_names.split()[-1]
                else:
                    parsed_data[s_name]["unclassified"] = {"count": int(reads), "percent": float(pct)}
            else:
                parsed_data[s_name]["assigned"][taxon_names] = {"count": int(reads), "percent": float(pct)}

        return taxo_rank, parsed_data

    def sum_sample_counts(self):
        """Sum counts across all samples for kaiju data"""

        # Sum the percentages for each taxa across all samples
        # Allows us to pick top-5 for each rank
        # Use percentages instead of counts so that deeply-sequences samples
        # are not unfairly over-represented
        first_pass_total = False

        for rank_name, data in self.kaiju_data.items():
            for s_name, samples_values in data.items():
                # perform sum at first level only
                if rank_name not in self.kaiju_total_pct:
                    self.kaiju_total_pct[rank_name] = dict()
                # If total did not already compute
                if s_name not in self.kaiju_sample_total_readcounts:
                    first_pass_total = True
                    self.kaiju_sample_total_readcounts[s_name] = 0

                if "assigned" in samples_values:
                    for classif, row in samples_values["assigned"].items():
                        if classif not in self.kaiju_total_pct[rank_name]:
                            self.kaiju_total_pct[rank_name][classif] = 0
                        self.kaiju_total_pct[rank_name][classif] += row["percent"]
                        if first_pass_total:
                            self.kaiju_sample_total_readcounts[s_name] += row["count"]
                    if first_pass_total:
                        if "cannot be assigned" in samples_values:
                            self.kaiju_sample_total_readcounts[s_name] += samples_values["cannot be assigned"]["count"]
                        if "unclassified" in samples_values:
                            self.kaiju_sample_total_readcounts[s_name] += samples_values["unclassified"]["count"]
                            self.kaiju_sample_unclassified[s_name] = samples_values["unclassified"]["count"]
                        first_pass_total = False

    def kaiju_stats_table(self):
        """Take the parsed stats from the Kaiju reports and add them to the
        basic stats table at the top of the report"""
        headers = {}
        taxo_ranks = sorted(self.kaiju_data)
        # print only assigned at phylum rank in general table or the first available.
        if len(taxo_ranks) >= 1 and "phylum" in taxo_ranks:
            general_taxo_rank = "Phylum"
        else:
            general_taxo_rank = taxo_ranks[0].capitalize()

        headers["% Assigned"] = {
            "title": "% Reads assigned",
            "description": f"Percentage of reads assigned at {general_taxo_rank} rank",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
        }
        headers["assigned"] = {
            "title": f"{config.read_count_prefix} Reads assigned",
            "description": "Number of reads assigned ({})  at {} rank".format(
                config.read_count_desc, general_taxo_rank
            ),
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "Blues",
        }
        headers["% Unclassified"] = {
            "title": "% Unclassified",
            "description": "Percentage of reads that were unclassified",
            "suffix": "%",
            "max": 100,
            "scale": "OrRd",
        }
        tdata = {}
        for s_name, d in self.kaiju_sample_total_readcounts.items():
            tdata[s_name] = {
                # Default values in case if data is not available for this sample at this rank
                "assigned": 0,
                "% Assigned": 0,
            }
            total = self.kaiju_sample_total_readcounts[s_name]
            if total == 0:
                continue
            tdata[s_name]["% Unclassified"] = self.kaiju_sample_unclassified[s_name] * 100 / total
            if s_name in self.kaiju_data[general_taxo_rank.lower()]:
                tdata[s_name]["assigned"] = (
                    total
                    - self.kaiju_sample_unclassified[s_name]
                    - self.kaiju_data[general_taxo_rank.lower()][s_name]["cannot be assigned"]["count"]
                )
                tdata[s_name]["% Assigned"] = tdata[s_name]["assigned"] * 100 / total
        self.general_stats_addcols(tdata, headers)

    def top_five_barplot(self):
        """Add a bar plot showing the top-5 from each taxa rank"""

        # ordered rank used
        rank_used = []
        for rank in self.taxonomic_ranks:
            if rank in self.kaiju_data:
                rank_used.append(rank)

        pconfig = {
            "id": "kaiju-topfive-plot",
            "title": "Kaiju: Top taxa",
            "ylab": "Number of reads",
            "cpswitch_counts_label": "Number of reads",
            "data_labels": [
                {"name": rank.capitalize(), "ylab": f"Number of Reads at {rank.capitalize()} rank"}
                for rank in rank_used
            ],
        }

        pd = []
        cats = []

        for rank_name in rank_used:
            rank_cats = {}
            rank_data = {}
            # Loop through the summed tax percentages to get the top 5 across all samples

            top_five_sorted_pct = sorted(self.kaiju_total_pct[rank_name].items(), key=lambda x: x[1], reverse=True)[0:5]

            # retrieve top five counts
            counts_shown = {}
            for classif, pct_sum in top_five_sorted_pct:
                rank_cats[classif] = {"name": classif}
                # Pull out counts for this rank + classif from each sample
                for s_name, d in self.kaiju_data[rank_name].items():
                    if s_name not in rank_data:
                        rank_data[s_name] = dict()
                    if s_name not in counts_shown:
                        counts_shown[s_name] = 0
                    if classif in d["assigned"]:
                        if classif not in rank_data[s_name]:
                            rank_data[s_name][classif] = 0
                        rank_data[s_name][classif] += d["assigned"][classif]["count"]
                        counts_shown[s_name] += d["assigned"][classif]["count"]

            # Add in unclassified/cannot assigned/other
            for s_name in self.all_s_name:
                if s_name in self.kaiju_data[rank_name]:
                    d = self.kaiju_data[rank_name][s_name]
                    rank_data[s_name]["not assigned"] = d["cannot be assigned"]["count"]
                    rank_data[s_name]["unclassified"] = d["unclassified"]["count"]
                    rank_data[s_name]["other"] = (
                        self.kaiju_sample_total_readcounts[s_name]
                        - counts_shown[s_name]
                        - d["cannot be assigned"]["count"]
                        - d["unclassified"]["count"]
                    )
                    # This should never happen... But it does sometimes if the total read count is a bit off
                    if rank_data[s_name]["other"] < 0:
                        log.debug(
                            "Found negative 'other' count for {} ({}): {}".format(
                                s_name, rank_name, rank_data[s_name]["other"]
                            )
                        )
                        rank_data[s_name]["other"] = 0

                else:
                    rank_data[s_name] = {}
                    rank_data[s_name]["other"] = 0
                    rank_data[s_name]["missing info"] = (
                        self.kaiju_sample_total_readcounts[s_name] - self.kaiju_sample_unclassified[s_name]
                    )
                    rank_data[s_name]["unclassified"] = self.kaiju_sample_unclassified[s_name]

            rank_cats["other"] = {"name": "Other", "color": "#a65628"}
            rank_cats["not assigned"] = {"name": "Cannot be assigned", "color": "#cccccc"}
            rank_cats["missing info"] = {"name": "Missing info", "color": "#979a9a"}
            rank_cats["unclassified"] = {"name": "Unclassified", "color": "#d4949c"}

            self.write_data_file(rank_data, "multiqc_kaiju_" + rank_name)
            cats.append(rank_cats)
            pd.append(rank_data)

        if len(cats) == 0:
            log.debug("No data for Kaiju top five barplot. Skipping.")
            return

        self.add_section(
            name="Top taxa",
            anchor="kaiju-topfive",
            description="The number of reads falling into the top 5 taxa across different ranks.",
            helptext="""
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The counts for these top five taxa are then plotted for each of the taxa ranks found in logs.
                The unclassified count is always shown across all taxa ranks.
                The 'Cannot be assigned' count correspond to reads classified but not at this taxa rank.

                The category _"Other"_ shows the difference between the above total assingned read count and the sum of the read counts
                in the top 5 taxa shown. This should cover all taxa _not_ in the top 5.
            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )
