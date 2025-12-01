import logging
import re

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The module supports outputs from sylphtax, that look like the following:

    ```tsv
    clade_name	relative_abundance	sequence_abundance	ANI (if strain-level)
    d__Bacteria	100.00010000000002	99.99999999999999	NA
    d__Bacteria|p__Bacillota	24.640800000000002	18.712699999999998	NA
    d__Bacteria|p__Bacillota_A	47.333499999999994	52.5969	NA
    ```

    A bar graph is generated that shows the relative abundance for each sample that
    fall into the top-10 categories for each taxa rank. The top categories are calculated
    by summing the relative abundances across all samples.

    The number of top categories to plot can be customized in the config file:

    ```yaml
    sylphtax:
      top_n: 10
    ```
    """

    def __init__(
        self,
        name="Sylph-tax",
        anchor="sylphtax",
        href=["https://sylph-docs.github.io/", "https://sylph-docs.github.io/sylph-tax/"],
        info="Taxonomic profiling of metagenomic reads.",
        doi="10.1038/s41587-024-02412-y",
    ):
        super(MultiqcModule, self).__init__(
            name=name,
            anchor=anchor,
            href=href,
            info=info,
            doi=doi,
        )
        # Taxonomic ranks: include Sylph’s domain/realm/strain if desired
        self.t_ranks = {
            "t": "Strain",
            "s": "Species",
            "g": "Genus",
            "f": "Family",
            "o": "Order",
            "c": "Class",
            "p": "Phylum",
            "k": "Kingdom",
            "d": "Domain",
            "r": "Realm",
            "u": "No Taxonomy",
        }

        self.top_n = getattr(config, "sylphtax", {}).get("top_n", 10)

        self.sylph_raw_data = dict()
        for f in self.find_log_files("sylphtax", filehandles=True):
            self.parse_logs(f)
            self.add_data_source(f)

        self.sylph_raw_data = self.ignore_samples(self.sylph_raw_data)

        if len(self.sylph_raw_data) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.sylph_raw_data)} reports")

        self.add_software_version(None)

        # Sum percentages across all samples, so that we can pick top species
        self.sylph_total_pct = dict()
        self.sum_sample_abundances()
        self.general_stats_cols()
        self.top_taxa_barplot()
        self.write_data_file(self.sylph_raw_data, f"multiqc_{self.anchor}")

    def parse_logs(self, f):
        """
        Parses a sylphtax .sylphmpa file to extract:
        1. A rank code, indicating:
            * (k)ingdom
            * (p)hylum
            * (c)lass
            * (o)rder
            * (f)amily
            * (g)enus
            * (s)pecies
            * s(t)rain
        2. The last clade name (taxonomy)
        3. Relative abundance of this clade (percentage)
        """
        regex_clade_rel = re.compile(
            r"^(?P<clade>[^\t#]+)\t(?P<rel_abundance>[+-]?(?:\d+(?:\.\d+)?|\.\d+)(?:[Ee][+-]?\d+)?)"
        )
        regex_last_tax_rank = re.compile(r"(?:^|[|;])([drkpcofgst])__([^|;]+)$")
        regex_comments = re.compile(r"^\s*#")
        regex_header = re.compile(r"^clade_name\t", re.IGNORECASE)

        data = []
        for line in f["f"]:
            if regex_comments.search(line) or regex_header.search(line):
                continue

            m = regex_clade_rel.match(line)
            if not m:
                log.debug(f"{f['s_name']}: Could not parse line: {line.strip()}")
                continue

            clade = m.group("clade")
            rel = float(m.group("rel_abundance"))

            # Handle bare NO_TAXONOMY rows (no rank code present)
            if clade == "NO_TAXONOMY":
                for rank in list("kcpofgsd"):
                    data.append(
                        {
                            "tax_rank": rank,
                            "taxonomy": "No Taxonomy (see strain)",
                            "rel_abundance": 100,
                        }
                    )
                continue

            # Extract last rank and taxonomy from clade path
            match_last_rank = regex_last_tax_rank.search(clade)
            if match_last_rank:
                row = {
                    "tax_rank": match_last_rank.group(1),
                    "taxonomy": match_last_rank.group(2),
                    "rel_abundance": rel,
                }
                data.append(row)
                continue

            # Fallback: single-token clade_name like "d__Bacteria"
            m_single = re.match(r"^([drkpcofgst])__([^\t|;]+)$", clade)
            if m_single:
                row = {
                    "tax_rank": m_single.group(1),
                    "taxonomy": m_single.group(2),
                    "rel_abundance": rel,
                }
                data.append(row)
            else:
                log.debug(f"{f['s_name']}: Could not parse clade_name: {clade}")

        self.sylph_raw_data[f["s_name"]] = data

    def sum_sample_abundances(self):
        """Sum relative abundance across all samples for sylph data"""

        # Sum the percentages for each taxa across all samples
        # Allows us to pick the top taxa for each rank
        for s_name, data in self.sylph_raw_data.items():
            for row in data:
                tax_rank = row["tax_rank"]
                taxonomy = row["taxonomy"]

                if tax_rank not in self.sylph_total_pct:
                    self.sylph_total_pct[tax_rank] = dict()

                if taxonomy not in self.sylph_total_pct[tax_rank]:
                    self.sylph_total_pct[tax_rank][taxonomy] = 0
                self.sylph_total_pct[tax_rank][taxonomy] += row["rel_abundance"]

    def general_stats_cols(self):
        # Find top taxa in the most specific non-strain rank available
        top_taxa = []
        top_rank_code = None
        top_rank_name = None

        for rank_code, rank_name in self.t_ranks.items():
            if rank_code == "t" or rank_code == "u":
                continue
            try:
                sorted_pct = sorted(
                    self.sylph_total_pct[rank_code].items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
                if sorted_pct:
                    top_taxa = [taxonomy for taxonomy, _ in sorted_pct[: self.top_n]]
                    top_rank_code = rank_code
                    top_rank_name = rank_name
                    break
            except KeyError:
                # Rank not present; try next
                pass

        # Fallbacks: strain first, then no taxonomy
        if not top_taxa:
            if "t" in self.sylph_total_pct:
                sorted_pct = sorted(
                    self.sylph_total_pct["t"].items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
                if sorted_pct:
                    top_taxa = [taxonomy for taxonomy, _ in sorted_pct[: self.top_n]]
                    top_rank_code = "t"
                    top_rank_name = self.t_ranks["t"]
            elif "u" in self.sylph_total_pct:
                sorted_pct = sorted(
                    self.sylph_total_pct["u"].items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
                if sorted_pct:
                    top_taxa = [taxonomy for taxonomy, _ in sorted_pct[: self.top_n]]
                    top_rank_code = "u"
                    top_rank_name = self.t_ranks["u"]

        # If still no taxa found, skip adding columns
        if not top_taxa or top_rank_code is None:
            log.warning("Sylph: No taxa found to populate General Stats")
            return

        # Column headers
        headers = dict()
        top_one_hkey = f"% {top_taxa[0]}"
        headers[top_one_hkey] = {
            "title": top_one_hkey,
            "description": "Percentage of reads that were the top {} over all samples - {}".format(
                top_rank_name, top_taxa[0]
            ),
            "suffix": "%",
            "max": 100,
            "scale": "PuBuGn",
        }
        headers["% Top"] = {
            "title": f"% Top {self.top_n} {top_rank_name}",
            "description": f"Percentage of reads that were classified by one of the top-{self.top_n} {top_rank_name} ({', '.join(top_taxa)})",
            "suffix": "%",
            "max": 100,
            "scale": "PuBu",
        }

        # Populate table data
        tdata = {}
        for s_name, d in self.sylph_raw_data.items():
            tdata[s_name] = {}
            for row in d:
                percent = row["rel_abundance"]
                if row["tax_rank"] == top_rank_code and row["taxonomy"] in top_taxa:
                    tdata[s_name]["% Top"] = percent + tdata[s_name].get("% Top", 0)
                if row["tax_rank"] == top_rank_code and row["taxonomy"] == top_taxa[0]:
                    tdata[s_name][top_one_hkey] = percent

            # Ensure presence of the single-top key even if absent in sample
            if top_one_hkey not in tdata[s_name]:
                tdata[s_name][top_one_hkey] = 0

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
                sorted_pct = sorted(
                    self.sylph_total_pct[rank_code].items(),
                    key=lambda x: x[1],
                    reverse=True,
                )
            except KeyError:
                # Taxa rank not found in this sample
                continue
            i = 0
            abundances_shown = {}
            for taxonomy, pct_sum in sorted_pct:
                # Add top n taxa for each rank
                i += 1
                if i > self.top_n:
                    # After top N, keep looping to sum up unclassified
                    continue
                rank_cats[taxonomy] = {"name": taxonomy}
                # Pull out abundances for this rank + classif from each sample
                for s_name, d in self.sylph_raw_data.items():
                    if s_name not in rank_data:
                        rank_data[s_name] = dict()
                    if s_name not in abundances_shown:
                        abundances_shown[s_name] = 0

                    for row in d:
                        if row["tax_rank"] == rank_code:
                            found_rank_codes.add(rank_code)
                            if row["taxonomy"] == taxonomy:
                                if taxonomy not in rank_data[s_name]:
                                    rank_data[s_name][taxonomy] = 0
                                rank_data[s_name][taxonomy] += row["rel_abundance"]
                                abundances_shown[s_name] += row["rel_abundance"]
                                rank_data[s_name]["other"] = 100 - abundances_shown[s_name]
            # Add in other - we presume from other species etc.
            for s_name, d in self.sylph_raw_data.items():
                # In case none of the top_n were in some sample:
                if s_name not in rank_data:
                    rank_data[s_name] = dict()
                if s_name not in abundances_shown:
                    abundances_shown[s_name] = 0
                rank_data[s_name]["other"] = 100 - abundances_shown[s_name]

                # This should never happen... But it does in Metaphlan at least if the total abundance is a bit off
                if rank_data[s_name]["other"] < 0:
                    log.debug(
                        "Found negative 'other' abundance for {} ({}): {}".format(
                            s_name, self.t_ranks[rank_code], rank_data[s_name]["other"]
                        )
                    )
                    rank_data[s_name]["other"] = 0
                # Quick fix to ensure that the "other" category is on end of bar plot
                rank_data[s_name]["zzz_other"] = rank_data[s_name].pop("other")

            rank_cats["zzz_other"] = {"name": "Other", "color": "#cccccc"}

            cats.append(rank_cats)
            pd.append(rank_data)

        pconfig = {
            "id": f"{self.anchor}-top-n-plot",
            "title": f"{self.name}: Top taxa",
            "ylab": "Relative Abundance",
            "data_labels": [v for k, v in self.t_ranks.items() if k in found_rank_codes],
            "cpswitch": False,
        }

        self.add_section(
            name="Top taxa",
            anchor=f"{self.anchor}-top-n",
            description=f"The relative abundance of reads falling into the top {self.top_n} taxa across different ranks.",
            helptext=f"""
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The relative abundance for these top {self.top_n} taxa are then plotted for each of the different taxa ranks.

                The category _"Other"_ shows the difference between 100% and the sum of the percent
                in the top {self.top_n} taxa shown. This should cover all taxa _not_ in the top {self.top_n}, +/- any rounding errors.
                Note that Sylph does not estimate the percent of unclassified reads, see [here](https://github.com/bluenote-1577/sylph/issues/49).
            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )
