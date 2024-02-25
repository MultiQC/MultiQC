""" MultiQC module to parse output from kraken """


import logging
import re

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """Kraken module"""

    def __init__(
        self,
        name="Kraken",
        anchor="kraken",
        href="https://ccb.jhu.edu/software/kraken/",
        info="is a taxonomic classification tool that uses exact k-mer matches to find the lowest common ancestor (LCA) of a given sequence.",
        doi="10.1186/gb-2014-15-3-r46",
        sp_key="kraken",
    ):
        super(MultiqcModule, self).__init__(
            name=name,
            anchor=anchor,
            href=href,
            info=info,
            doi=doi,
        )
        self.t_ranks = {
            "S": "Species",
            "G": "Genus",
            "F": "Family",
            "O": "Order",
            "C": "Class",
            "P": "Phylum",
            "K": "Kingdom",
            "D": "Domain",
            "R": "Root",
            "U": "Unclassified",
        }

        self.top_n = getattr(config, "kraken", {}).get("top_n", 5)

        # Find and load any kraken reports
        self.kraken_raw_data = dict()
        new_report_present = False
        for f in self.find_log_files(sp_key, filehandles=True):
            log_is_new = self.log_is_new(f)
            f["f"].seek(0)
            if not log_is_new:
                self.parse_logs(f)
            else:
                new_report_present = True
                self.parse_logs_minimizer(f)
            self.add_data_source(f)

        self.kraken_raw_data = self.ignore_samples(self.kraken_raw_data)

        if len(self.kraken_raw_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.kraken_raw_data)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(self.kraken_raw_data, f"multiqc_{self.anchor}")

        # Sum counts across all samples, so that we can pick top species
        self.kraken_total_pct = dict()
        self.kraken_total_counts = dict()
        self.kraken_sample_total_readcounts = dict()
        self.sample_total_readcounts()
        self.sum_sample_counts()

        self.general_stats_cols()
        self.top_taxa_barplot()
        if new_report_present:
            self.top_taxa_duplication_heatmap()

    @staticmethod
    def log_is_new(f):
        """Check which version of Kraken report file is used

        If 6 fields are used, it's the 'old' log (without distinct minimizer)
        if 8 fields, the new log experimental log (with distinct minimizer)
        """

        for line in f["f"]:
            if len(line.split()) > 6:
                return True
            else:
                return False
        return False

    def parse_logs(self, f):
        """
        Parses a kraken report output file

        1. Percentage of fragments covered by the clade rooted at this taxon
        2. Number of fragments covered by the clade rooted at this taxon
        3. Number of fragments assigned directly to this taxon
        4. A rank code, indicating:
            * (U)nclassified
            * (R)oot
            * (D)omain
            * (K)ingdom
            * (P)hylum
            * (C)lass
            * (O)rder
            * (F)amily
            * (G)enus
            * (S)pecies
           Taxa that are not at any of these 10 ranks have a rank code that is
           formed by using the rank code of the closest ancestor rank with
           a number indicating the distance from that rank.  E.g., "G2" is a
           rank code indicating a taxon is between genus and species and the
           grandparent taxon is at the genus rank.
        5. NCBI taxonomic ID number
        6. Indented scientific name
        """

        # Search regexes for stats
        k2_regex = re.compile(r"^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t([\dUDKRPCOFGS-]{1,3})\t(\d+)(\s+)(.+)")
        data = []
        for line in f["f"]:
            match = k2_regex.search(line)
            if match:
                row = {
                    "percent": float(match.group(1)),
                    "counts_rooted": int(match.group(2)),
                    "counts_direct": int(match.group(3)),
                    "rank_code": match.group(4),
                    "tax_id": int(match.group(5)),
                    "num_spaces": len(match.group(6)),
                    "classif": match.group(7),
                }
                data.append(row)

        self.kraken_raw_data[f["s_name"]] = data

    def parse_logs_minimizer(self, f):
        """
        Parses a kraken report output file

        1. Percentage of fragments covered by the clade rooted at this taxon
        2. Number of fragments covered by the clade rooted at this taxon
        3. Number of fragments assigned directly to this taxon
        4. Number of minimizers in read data associated with this taxon (new)
        5. An estimate of the number of distinct minimizers in read data
           associated with this taxon (new)
        6. A rank code, indicating:
            * (U)nclassified
            * (R)oot
            * (D)omain
            * (K)ingdom
            * (P)hylum
            * (C)lass
            * (O)rder
            * (F)amily
            * (G)enus
            * (S)pecies
           Taxa that are not at any of these 10 ranks have a rank code that is
           formed by using the rank code of the closest ancestor rank with
           a number indicating the distance from that rank.  E.g., "G2" is a
           rank code indicating a taxon is between genus and species and the
           grandparent taxon is at the genus rank.
        7. NCBI taxonomic ID number
        8. Indented scientific name
        """

        def duplication(total, distinct):
            """Get duplication factor

            Args:
                total (int): Total number of elements
                distinct (int): Number of distinct elements
            """
            try:
                res = float(total) / distinct
            except ZeroDivisionError:
                res = 0.0
            return res

        # Search regexes for stats
        k2_regex = re.compile(
            r"^\s{0,2}(\d{1,3}\.\d{1,2})\t(\d+)\t(\d+)\t(\d+)\t(\d+)\t([URDKPCOFGS-]\d{0,2})\t(\d+)(\s+)(.+)"
        )
        data = []
        for line in f["f"]:
            match = k2_regex.search(line)
            if match:
                row = {
                    "percent": float(match.group(1)),
                    "counts_rooted": int(match.group(2)),
                    "counts_direct": int(match.group(3)),
                    "minimizer": int(match.group(4)),
                    "minimizer_distinct": int(match.group(5)),
                    "minimizer_duplication": duplication(int(match.group(4)), int(match.group(5))),
                    "rank_code": match.group(6),
                    "tax_id": int(match.group(7)),
                    "num_spaces": len(match.group(8)),
                    "classif": match.group(9),
                }
                data.append(row)
            else:
                log.debug(f"{f['s_name']}: Could not parse line: {line}")

        self.kraken_raw_data[f["s_name"]] = data

    def sample_total_readcounts(self):
        """Compute the total read counts for each sample"""

        total_all_samples = 0
        for s_name, data in self.kraken_raw_data.items():
            self.kraken_sample_total_readcounts[s_name] = 0
            for row in data:
                self.kraken_sample_total_readcounts[s_name] += row["counts_direct"]
            total_all_samples += self.kraken_sample_total_readcounts[s_name]

        # Check that we had some counts for some samples, exit if not
        if total_all_samples == 0:
            log.warning("No samples had any reads")
            raise ModuleNoSamplesFound

    def sum_sample_counts(self):
        """Sum counts across all samples for kraken data"""

        # Sum the percentages for each taxa across all samples
        # Allows us to pick the top taxa for each rank
        # Use percentages instead of counts so that deeply-sequences samples
        # are not unfairly over-represented
        for s_name, data in self.kraken_raw_data.items():
            for row in data:
                # Convenience vars that are easier to read
                rank_code = row["rank_code"]
                classif = row["classif"]

                # Skip anything that doesn't exactly fit a tax rank level
                if row["rank_code"] == "-" or any(c.isdigit() for c in row["rank_code"]):
                    continue

                if rank_code not in self.kraken_total_pct:
                    self.kraken_total_pct[rank_code] = dict()
                    self.kraken_total_counts[rank_code] = dict()

                if classif not in self.kraken_total_pct[rank_code]:
                    self.kraken_total_pct[rank_code][classif] = 0
                    self.kraken_total_counts[rank_code][classif] = 0
                try:
                    self.kraken_total_pct[rank_code][classif] += (
                        row["counts_rooted"] / self.kraken_sample_total_readcounts[s_name]
                    )
                except ZeroDivisionError:
                    pass
                self.kraken_total_counts[rank_code][classif] += row["counts_rooted"]

    def general_stats_cols(self):
        """Add a couple of columns to the General Statistics table"""

        # Get top taxa in most specific taxa rank that we have
        top_taxa = []
        top_rank_code = None
        top_rank_name = None
        for rank_code, rank_name in self.t_ranks.items():
            try:
                sorted_pct = sorted(self.kraken_total_pct[rank_code].items(), key=lambda x: x[1], reverse=True)
                for classif, pct_sum in sorted_pct[: self.top_n]:
                    top_taxa.append(classif)
                top_rank_code = rank_code
                top_rank_name = rank_name
                break
            except KeyError:
                # No species-level data found etc
                pass

        # Column headers
        headers = dict()

        top_one = None

        # don't include top-N % in general stats if all is unclassified.
        # unclassified is included separately, so also don't include twice
        if top_rank_code != "U":
            top_one = f"% {top_taxa[0]}"
            headers["pct_top_one"] = {
                "title": top_one,
                "description": "Percentage of reads that were the top {} over all samples - {}".format(
                    top_rank_name, top_taxa[0]
                ),
                "suffix": "%",
                "max": 100,
                "scale": "PuBuGn",
            }
            headers["pct_top_n"] = {
                "title": f"% Top {self.top_n} {top_rank_name}",
                "description": f"Percentage of reads that were classified by one of the top-{self.top_n} {top_rank_name} ({', '.join(top_taxa)})",
                "suffix": "%",
                "max": 100,
                "scale": "PuBu",
            }

        headers["pct_unclassified"] = {
            "title": "% Unclassified",
            "description": "Percentage of reads that were unclassified",
            "suffix": "%",
            "max": 100,
            "scale": "OrRd",
        }

        # Get table data
        tdata = {}
        for s_name, d in self.kraken_raw_data.items():
            tdata[s_name] = {}
            for row in d:
                try:
                    percent = (row["counts_rooted"] / self.kraken_sample_total_readcounts[s_name]) * 100.0
                except ZeroDivisionError:
                    percent = 0
                if row["rank_code"] == "U":
                    tdata[s_name]["pct_unclassified"] = percent
                if row["rank_code"] == top_rank_code and row["classif"] in top_taxa:
                    tdata[s_name]["pct_top_n"] = percent + tdata[s_name].get("pct_top_n", 0)
                if row["rank_code"] == top_rank_code and row["classif"] == top_taxa[0]:
                    tdata[s_name]["pct_top_one"] = percent

            if top_one is not None and "pct_top_one" not in tdata[s_name]:
                tdata[s_name]["pct_top_one"] = 0

        self.general_stats_addcols(tdata, headers)

    def top_taxa_barplot(self):
        """Add a bar plot showing the top-N from each taxa rank"""

        pd = []
        cats = list()
        # Keeping track of encountered codes to display only tabs with available data
        found_rank_codes = set()

        for rank_code in self.t_ranks:
            rank_cats = dict()
            rank_data = dict()

            # Loop through the summed tax percentages to get the top-N across all samples
            try:
                sorted_pct = sorted(self.kraken_total_pct[rank_code].items(), key=lambda x: x[1], reverse=True)
            except KeyError:
                # Taxa rank not found in this sample
                continue
            i = 0
            counts_shown = {}
            for classif, pct_sum in sorted_pct:
                i += 1
                if i > self.top_n:
                    break
                rank_cats[classif] = {"name": classif}
                # Pull out counts for this rank + classif from each sample
                for s_name, d in self.kraken_raw_data.items():
                    if s_name not in rank_data:
                        rank_data[s_name] = dict()
                    if s_name not in counts_shown:
                        counts_shown[s_name] = 0

                    for row in d:
                        if row["rank_code"] == rank_code:
                            found_rank_codes.add(rank_code)
                            # unclassified are handled separately
                            if row["rank_code"] != "U":
                                if row["classif"] == classif:
                                    if classif not in rank_data[s_name]:
                                        rank_data[s_name][classif] = 0
                                    rank_data[s_name][classif] += row["counts_rooted"]
                                    counts_shown[s_name] += row["counts_rooted"]

            # Add in unclassified reads and "other" - we presume from other species etc.
            for s_name, d in self.kraken_raw_data.items():
                for row in d:
                    if row["rank_code"] == "U":
                        rank_data[s_name]["U"] = row["counts_rooted"]
                        counts_shown[s_name] += row["counts_rooted"]
                rank_data[s_name]["other"] = self.kraken_sample_total_readcounts[s_name] - counts_shown[s_name]

                # This should never happen... But it does sometimes if the total read count is a bit off
                if rank_data[s_name]["other"] < 0:
                    log.debug(
                        "Found negative 'other' count for {} ({}): {}".format(
                            s_name, self.t_ranks[rank_code], rank_data[s_name]["other"]
                        )
                    )
                    rank_data[s_name]["other"] = 0

            rank_cats["other"] = {"name": "Other", "color": "#cccccc"}
            rank_cats["U"] = {"name": "Unclassified", "color": "#d4949c"}

            cats.append(rank_cats)
            pd.append(rank_data)

        pconfig = {
            "id": f"{self.anchor}-top-n-plot",
            "title": f"{self.name}: Top taxa",
            "ylab": "Number of fragments",
            "data_labels": [v for k, v in self.t_ranks.items() if k in found_rank_codes],
        }

        self.add_section(
            name="Top taxa",
            anchor=f"{self.anchor}-top-n",
            description=f"The number of reads falling into the top {self.top_n} taxa across different ranks.",
            helptext=f"""
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The counts for these top {self.top_n} taxa are then plotted for each of the 9 different taxa ranks.
                The unclassified count is always shown across all taxa ranks.

                The total number of reads is approximated by dividing the number of `unclassified` reads by the percentage of
                the library that they account for.
                Note that this is only an approximation, and that kraken percentages don't always add to exactly 100%.

                The category _"Other"_ shows the difference between the above total read count and the sum of the read counts
                in the top {self.top_n} taxa shown + unclassified. This should cover all taxa _not_ in the top {self.top_n}, +/- any rounding errors.

                Note that any taxon that does not exactly fit a taxon rank (eg. `-` or `G2`) is ignored.
            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )

    def top_taxa_duplication_heatmap(self):
        """Add a heatmap showing the minimizer duplication of the top taxa"""

        duplication = list()
        pconfig = {
            "id": f"{self.anchor}-top-duplication_plot",
            "title": f"{self.name}: Top {self.top_n} species duplication",
            "square": False,
            "xcats_samples": False,
            "angled_xticks": False,
        }

        rank_code = "S"
        rank_data = dict()
        # Loop through the summed tax percentages to get the top taxa across all samples
        try:
            sorted_pct = sorted(self.kraken_total_pct[rank_code].items(), key=lambda x: x[1], reverse=True)
        except KeyError:
            log.debug("Taxa rank not found, skipping Taxa duplication heatmap")
            return

        i = 0
        counts_shown = {}
        showed_warning = False
        for classif, pct_sum in sorted_pct:
            i += 1
            if i > self.top_n:
                break
            # Pull out counts for this rank + classif from each sample
            for s_name, d in self.kraken_raw_data.items():
                if s_name not in rank_data:
                    rank_data[s_name] = dict()
                if s_name not in counts_shown:
                    counts_shown[s_name] = 0

                if classif not in rank_data[s_name]:
                    rank_data[s_name][classif] = None

                try:
                    row = next(row for row in d if row["rank_code"] == rank_code and row["classif"] == classif)
                except StopIteration:
                    # if nothing is found at the rank + classification, leave as 0
                    continue

                try:
                    rank_data[s_name][classif] = row["minimizer_duplication"]
                except KeyError:
                    del rank_data[s_name]
                    if not showed_warning:
                        log.warning("Kraken2 reports of different versions were found")
                        showed_warning = True

        # Strip empty samples
        for sample, vals in dict(rank_data).items():
            if len(vals) == 0:
                del rank_data[sample]

        # Build data structures for heatmap
        ylabels = list(rank_data.keys())
        xlabels = list(rank_data[ylabels[0]].keys())
        for sample in rank_data:
            duplication.append(list(rank_data[sample].values()))

        self.add_section(
            name="Duplication rate of top species",
            anchor=f"{self.anchor}-duplication-topfive",
            description=f"The duplication rate of minimizer falling into the top {self.top_n} species",
            helptext=f"""
                To make this plot, the minimizer duplication rate is computed for the top {self.top_n} most abundant species in all samples.

                The minimizer duplication rate is defined as: `duplication rate = (total number of minimizers / number of distinct minimizers)`

                A low coverage and high duplication rate (`>> 1`) is often sign of read stacking, which probably indicates of false positive hit.
            """,
            plot=heatmap.plot(duplication, xlabels, ylabels, pconfig),
        )
