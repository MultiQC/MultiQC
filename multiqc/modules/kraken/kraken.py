import logging
import re
from collections import defaultdict
from typing import Dict, Union, List

from multiqc import config
from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, heatmap

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """
    The MultiQC module supports outputs from both Kraken and Kraken 2.

    It works with report files generated using the `--report` flag, that look like the following:

    ```ts
    11.66	98148	98148	U	0	unclassified
    88.34	743870	996	-	1	root
    88.22	742867	0	-	131567	  cellular organisms
    88.22	742866	2071	D	2	    Bacteria
    87.95	740514	2914	P	1239	      Firmicutes
    ```

    A bar graph is generated that shows the number of fragments for each sample that
    fall into the top-5 categories for each taxa rank. The top categories are calculated
    by summing the library percentages across all samples.

    The number of top categories to plot can be customized in the config file:

    ```yaml
    kraken:
      top_n: 5
    ```
    """

    T_RANKS = {
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

    TOP_N = getattr(config, "kraken", {}).get("top_n", 5)

    def __init__(
        self,
        name="Kraken",
        anchor="kraken",
        href="https://ccb.jhu.edu/software/kraken/",
        info="Taxonomic classification tool that uses exact k-mer matches to find the lowest common ancestor "
        "(LCA) of a given sequence.",
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

        # Find and load any kraken reports
        raw_rows_by_sample: Dict[str, List[Dict[str, Union[int, str, float]]]] = dict()
        new_report_present = False

        for f in self.find_log_files(sp_key, filehandles=True):
            if not log_is_new(f):
                raw_rows_by_sample[f["s_name"]] = parse_logs(f)
            else:
                new_report_present = True
                raw_rows_by_sample[f["s_name"]] = parse_logs_minimizer(f)
            self.add_data_source(f)

        raw_rows_by_sample = self.ignore_samples(raw_rows_by_sample)
        if len(raw_rows_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"{name + ': ' if name != 'Kraken' else ''}Found {len(raw_rows_by_sample)} reports")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(raw_rows_by_sample, f"multiqc_{self.anchor}")

        # Sum counts across all samples, so that we can pick top species
        total_cnt_by_sample = self.sample_total_readcounts(raw_rows_by_sample)
        cnt_by_classif_by_rank, pct_by_classif_by_rank = sum_sample_counts(raw_rows_by_sample, total_cnt_by_sample)

        self.general_stats_cols(raw_rows_by_sample, total_cnt_by_sample, pct_by_classif_by_rank)
        self.top_taxa_barplot(raw_rows_by_sample, total_cnt_by_sample, pct_by_classif_by_rank)
        if new_report_present:
            self.top_taxa_duplication_heatmap(raw_rows_by_sample, pct_by_classif_by_rank)

    def sample_total_readcounts(
        self, rows_by_sample: Dict[str, List[Dict[str, Union[str, int, float]]]]
    ) -> Dict[str, int]:
        """Compute the total read counts for each sample"""
        total_cnt_by_sample: Dict[str, int] = dict()

        _total_all_samples = 0
        # Take the unassigned counts (line 1) and counts assigned to root (line 2) for each sample
        for s_name, row in rows_by_sample.items():
            # The 2nd column here just contains the number of unassigned reads
            unassigned_counts = int(row[0]["counts_rooted"])
            # the 2nd column in other rows contains the number of reads mapped to this taxa
            assigned_counts = int(row[1]["counts_rooted"]) if len(row) > 1 else 0
            total_cnt_by_sample[s_name] = unassigned_counts + assigned_counts
            _total_all_samples += total_cnt_by_sample[s_name]

        # Check that we had some counts for some samples, exit if not
        if _total_all_samples == 0:
            log.warning("No samples had any reads")
            raise ModuleNoSamplesFound

        return total_cnt_by_sample

    def general_stats_cols(
        self,
        rows_by_sample: Dict[str, List[Dict[str, Union[str, int, float]]]],
        total_cnt_by_sample: Dict[str, int],
        total_pct_by_rank_by_classif: Dict[str, Dict[str, float]],
    ):
        """Add a couple of columns to the General Statistics table"""

        # Get top taxa in most specific taxa rank that we have
        top_taxa = []
        top_rank_code = None
        top_rank_name = None
        for rank_code, rank_name in MultiqcModule.T_RANKS.items():
            if rank_code in total_pct_by_rank_by_classif:
                sorted_pct = sorted(total_pct_by_rank_by_classif[rank_code].items(), key=lambda x: x[1], reverse=True)
                for classif, pct_sum in sorted_pct[: MultiqcModule.TOP_N]:
                    top_taxa.append(classif)
                top_rank_code = rank_code
                top_rank_name = rank_name
                break

        if not top_taxa:
            log.error("No taxa found")
            return

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
                "title": f"% Top {MultiqcModule.TOP_N} {top_rank_name}",
                "description": f"Percentage of reads that were classified by one of the top-{MultiqcModule.TOP_N} {top_rank_name} ({', '.join(top_taxa)})",
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
        table_data_by_sample: Dict[str, Dict[str, float]] = {}
        for s_name, d in rows_by_sample.items():
            table_data_by_sample[s_name] = {}
            for row in d:
                if total_cnt_by_sample[s_name] != 0:
                    percent = (int(row["counts_rooted"]) / total_cnt_by_sample[s_name]) * 100.0
                else:
                    percent = 0
                if row["rank_code"] == "U":
                    table_data_by_sample[s_name]["pct_unclassified"] = percent
                if row["rank_code"] == top_rank_code and row["classif"] in top_taxa:
                    table_data_by_sample[s_name]["pct_top_n"] = percent + table_data_by_sample[s_name].get(
                        "pct_top_n", 0
                    )
                if row["rank_code"] == top_rank_code and row["classif"] == top_taxa[0]:
                    table_data_by_sample[s_name]["pct_top_one"] = percent

            if top_one is not None and "pct_top_one" not in table_data_by_sample[s_name]:
                table_data_by_sample[s_name]["pct_top_one"] = 0

        self.general_stats_addcols(table_data_by_sample, headers)

    def top_taxa_barplot(
        self,
        rows_by_sample: Dict[str, List[Dict[str, Union[str, int, float]]]],
        total_cnt_by_sample: Dict[str, int],
        total_pct_by_rank_by_classif: Dict[str, Dict[str, float]],
    ):
        """Add a bar plot showing the top-N from each taxa rank"""

        pd = []
        cats = list()
        # Keeping track of encountered codes to display only tabs with available data
        found_rank_codes = set()

        for rank_code in MultiqcModule.T_RANKS:
            rank_cats = dict()
            data_by_classif_by_sample: Dict[str, Dict[str, int]] = dict()

            # Loop through the summed tax percentages to get the top-N across all samples
            if rank_code not in total_pct_by_rank_by_classif:
                # Taxa rank not found in this sample
                continue

            sorted_pct = sorted(total_pct_by_rank_by_classif[rank_code].items(), key=lambda x: x[1], reverse=True)
            i = 0
            counts_shown: Dict[str, int] = dict()
            for classif, pct_sum in sorted_pct:
                i += 1
                if i > MultiqcModule.TOP_N:
                    break
                rank_cats[classif] = {"name": classif}
                # Pull out counts for this rank + classif from each sample
                for s_name, rows in rows_by_sample.items():
                    if s_name not in data_by_classif_by_sample:
                        data_by_classif_by_sample[s_name] = dict()
                    if s_name not in counts_shown:
                        counts_shown[s_name] = 0

                    for row in rows:
                        if row["rank_code"] == rank_code:
                            found_rank_codes.add(rank_code)
                            # unclassified are handled separately
                            if row["rank_code"] != "U":
                                if row["classif"] == classif:
                                    if classif not in data_by_classif_by_sample[s_name]:
                                        data_by_classif_by_sample[s_name][classif] = 0
                                    data_by_classif_by_sample[s_name][classif] += int(row["counts_rooted"])
                                    counts_shown[s_name] += int(row["counts_rooted"])

            # Add in unclassified reads and "other" - we presume from other species etc.
            for s_name, rows in rows_by_sample.items():
                for row in rows:
                    if row["rank_code"] == "U":
                        data_by_classif_by_sample[s_name]["U"] = int(row["counts_rooted"])
                        counts_shown[s_name] += int(row["counts_rooted"])
                data_by_classif_by_sample[s_name]["other"] = total_cnt_by_sample[s_name] - counts_shown[s_name]

                # This should never happen... But it does sometimes if the total read count is a bit off
                if data_by_classif_by_sample[s_name]["other"] < 0:
                    log.debug(
                        "Found negative 'other' count for {} ({}): {}".format(
                            s_name, MultiqcModule.T_RANKS[rank_code], data_by_classif_by_sample[s_name]["other"]
                        )
                    )
                    data_by_classif_by_sample[s_name]["other"] = 0

            rank_cats["other"] = {"name": "Other", "color": "#cccccc"}
            rank_cats["U"] = {"name": "Unclassified", "color": "#d4949c"}

            cats.append(rank_cats)
            pd.append(data_by_classif_by_sample)

        pconfig = {
            "id": f"{self.anchor}-top-n-plot",
            "title": f"{self.name}: Top taxa",
            "ylab": "Number of fragments",
            "data_labels": [v for k, v in MultiqcModule.T_RANKS.items() if k in found_rank_codes],
        }

        self.add_section(
            name="Top taxa",
            anchor=f"{self.anchor}-top-n",
            description=f"The number of reads falling into the top {MultiqcModule.TOP_N} taxa across different ranks.",
            helptext=f"""
                To make this plot, the percentage of each sample assigned to a given taxa is summed across all samples.
                The counts for these top {MultiqcModule.TOP_N} taxa are then plotted for each of the 9 different taxa ranks.
                The unclassified count is always shown across all taxa ranks.

                The total number of reads is approximated by dividing the number of `unclassified` reads by the percentage of
                the library that they account for.
                Note that this is only an approximation, and that kraken percentages don't always add to exactly 100%.

                The category _"Other"_ shows the difference between the above total read count and the sum of the read counts
                in the top {MultiqcModule.TOP_N} taxa shown + unclassified. This should cover all taxa _not_ in the top {MultiqcModule.TOP_N}, +/- any rounding errors.

                Note that any taxon that does not exactly fit a taxon rank (eg. `-` or `G2`) is ignored.
            """,
            plot=bargraph.plot(pd, cats, pconfig),
        )

    def top_taxa_duplication_heatmap(
        self,
        rows_by_sample: Dict[str, List[Dict[str, Union[str, int, float]]]],
        total_pct_by_rank_by_classif: Dict[str, Dict[str, float]],
    ):
        """Add a heatmap showing the minimizer duplication of the top taxa"""

        duplication = list()
        pconfig = {
            "id": f"{self.anchor}-top-duplication_plot",
            "title": f"{self.name}: Top {MultiqcModule.TOP_N} species duplication",
            "square": False,
            "xcats_samples": False,
            "angled_xticks": False,
        }

        rank_code = "S"
        rank_data_by_classif_by_sample: Dict[str, Dict[str, Union[int, None]]] = dict()
        # Loop through the summed tax percentages to get the top taxa across all samples
        if rank_code not in total_pct_by_rank_by_classif:
            log.debug(f"Taxa rank {rank_code} not found, skipping taxa duplication heatmap")
            return

        sorted_pct = sorted(total_pct_by_rank_by_classif[rank_code].items(), key=lambda x: x[1], reverse=True)
        i = 0
        counts_shown = {}
        showed_warning = False
        for classif, pct_sum in sorted_pct:
            i += 1
            if i > MultiqcModule.TOP_N:
                break
            # Pull out counts for this rank + classif from each sample
            for s_name, rows in rows_by_sample.items():
                if s_name not in rank_data_by_classif_by_sample:
                    rank_data_by_classif_by_sample[s_name] = dict()
                if s_name not in counts_shown:
                    counts_shown[s_name] = 0

                if classif not in rank_data_by_classif_by_sample[s_name]:
                    rank_data_by_classif_by_sample[s_name][classif] = None

                try:
                    row = next(row for row in rows if row["rank_code"] == rank_code and row["classif"] == classif)
                except StopIteration:
                    # if nothing is found at the rank + classification, leave as 0
                    continue

                try:
                    rank_data_by_classif_by_sample[s_name][classif] = int(row["minimizer_duplication"])
                except KeyError:
                    del rank_data_by_classif_by_sample[s_name]
                    if not showed_warning:
                        log.warning("Kraken2 reports of different versions were found")
                        showed_warning = True

        # Strip empty samples
        for sample, vals in dict(rank_data_by_classif_by_sample).items():
            if len(vals) == 0:
                del rank_data_by_classif_by_sample[sample]

        # Build data structures for heatmap
        ylabels = list(rank_data_by_classif_by_sample.keys())
        xlabels = list(rank_data_by_classif_by_sample[ylabels[0]].keys())
        for sample in rank_data_by_classif_by_sample:
            duplication.append(list(rank_data_by_classif_by_sample[sample].values()))

        self.add_section(
            name="Duplication rate of top species",
            anchor=f"{self.anchor}-duplication-topfive",
            description=f"The duplication rate of minimizer falling into the top {MultiqcModule.TOP_N} species",
            helptext=f"""
                To make this plot, the minimizer duplication rate is computed for the top {MultiqcModule.TOP_N} most abundant species in all samples.

                The minimizer duplication rate is defined as: `duplication rate = (total number of minimizers / number of distinct minimizers)`

                A low coverage and high duplication rate (`>> 1`) is often sign of read stacking, which probably indicates of false positive hit.
            """,
            plot=heatmap.plot(duplication, xlabels, ylabels, pconfig),
        )


def log_is_new(f) -> bool:
    """Check which version of Kraken report file is used

    If 6 fields are used, it's the 'old' log (without distinct minimizer)
    if 8 fields, the new log experimental log (with distinct minimizer)
    """

    result = False
    for line in f["f"]:
        if len(line.split()) > 6:
            result = True
            break
    f["f"].seek(0)
    return result


def parse_logs(f) -> List[Dict[str, Union[str, float, int]]]:
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

    return data


def parse_logs_minimizer(f) -> List[Dict[str, Union[str, int, float]]]:
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

    def duplication(total: int, distinct: int) -> float:
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
    data: List[Dict[str, Union[str, int, float]]] = []
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

    return data


def sum_sample_counts(
    rows_by_sample: Dict[str, List[Dict[str, Union[str, int, float]]]], total_cnt_by_sample: Dict[str, int]
):
    """Sum counts across all samples for kraken data"""

    cnt_by_classif_by_rank: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
    pct_by_classif_by_rank: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))

    # Sum the percentages for each taxon across all samples
    # Allows us to pick the top taxa for each rank
    # Use percentages instead of counts so that deeply-sequences samples
    # are not unfairly over-represented
    for s_name, rows in rows_by_sample.items():
        for row in rows:
            # Convenience vars that are easier to read
            rank_code = str(row["rank_code"])
            classif = str(row["classif"])

            # Skip anything that doesn't exactly fit a tax rank level
            if rank_code == "-" or any(c.isdigit() for c in rank_code):
                continue

            cnt_by_classif_by_rank[rank_code][classif] += int(row["counts_rooted"])
            if total_cnt_by_sample[s_name] != 0:
                pct_by_classif_by_rank[rank_code][classif] += int(row["counts_rooted"]) / total_cnt_by_sample[s_name]
    return cnt_by_classif_by_rank, pct_by_classif_by_rank
