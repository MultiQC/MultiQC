import logging
from collections import defaultdict
from typing import Dict, Union, List, Tuple

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

        total_cnt_by_sample: Dict[str, int] = dict()
        cnt_by_top_taxon_by_rank_by_sample: Dict[str, Dict[str, Dict[str, int]]] = dict()
        species_minimizer_dup_by_top_taxon_by_sample: Dict[str, Dict[str, float]] = dict()

        for f in self.find_log_files(sp_key, filehandles=True):
            sample_cnt_by_taxon_by_rank, min_dup_by_by_rank = parse_logs(f)

            # Sum the unassigned counts (line 1) and counts assigned to root (line 2) for each sample
            total_cnt = (
                sample_cnt_by_taxon_by_rank.get("U", {}).get("unclassified", 0)  # bracken doesn't have unclassified
                + sample_cnt_by_taxon_by_rank["R"]["root"]
            )
            if total_cnt == 0:
                log.warning(f"No reads found in {f['fn']}")
                continue

            self.add_data_source(f)
            if f["s_name"] in total_cnt_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

            total_cnt_by_sample[f["s_name"]] = total_cnt
            cnt_by_top_taxon_by_rank_by_sample[f["s_name"]] = sample_cnt_by_taxon_by_rank
            if min_dup_by_by_rank:
                species_minimizer_dup_by_top_taxon_by_sample[f["s_name"]] = min_dup_by_by_rank

        total_cnt_by_sample = self.ignore_samples(total_cnt_by_sample)
        if len(total_cnt_by_sample) == 0:
            raise ModuleNoSamplesFound
        log.info(f"{name + ': ' if name != 'Kraken' else ''}Found {len(total_cnt_by_sample)} reports")

        cnt_by_top_taxon_by_rank_by_sample = self.ignore_samples(cnt_by_top_taxon_by_rank_by_sample)
        species_minimizer_dup_by_top_taxon_by_sample = self.ignore_samples(species_minimizer_dup_by_top_taxon_by_sample)

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        self.write_data_file(cnt_by_top_taxon_by_rank_by_sample, f"multiqc_{self.anchor}")

        pct_by_top_taxon_by_rank: Dict[str, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for s_name, cnt_by_taxon_by_rank in cnt_by_top_taxon_by_rank_by_sample.items():
            for rank_code, cnt_by_taxon in cnt_by_taxon_by_rank.items():
                for taxon, count in cnt_by_taxon.items():
                    pct_by_top_taxon_by_rank[rank_code][taxon] += count / total_cnt_by_sample[s_name]

        self.general_stats_cols(total_cnt_by_sample, pct_by_top_taxon_by_rank, cnt_by_top_taxon_by_rank_by_sample)
        self.top_taxa_barplot(total_cnt_by_sample, pct_by_top_taxon_by_rank, cnt_by_top_taxon_by_rank_by_sample)
        if species_minimizer_dup_by_top_taxon_by_sample:
            self.top_taxa_duplication_heatmap(pct_by_top_taxon_by_rank, species_minimizer_dup_by_top_taxon_by_sample)

    def sample_total_readcounts(
        self, rows_by_sample: Dict[str, List[Dict[str, Union[str, int, float]]]]
    ) -> Dict[str, int]:
        """Compute the total read counts for each sample"""
        total_cnt_by_sample: Dict[str, int] = dict()

        _total_all_samples = 0
        # Take the unassigned counts (line 1) and counts assigned to root (line 2) for each sample
        for s_name, rows in rows_by_sample.items():
            # The 2nd column here just contains the number of unassigned reads
            unassigned_counts = int(rows[0]["counts_rooted"])
            # the 2nd column in other rows contains the number of reads mapped to this taxa
            assigned_counts = int(rows[1]["counts_rooted"]) if len(rows) > 1 else 0
            total_cnt_by_sample[s_name] = unassigned_counts + assigned_counts
            _total_all_samples += total_cnt_by_sample[s_name]

        # Check that we had some counts for some samples, exit if not
        if _total_all_samples == 0:
            log.warning("No samples had any reads")
            raise ModuleNoSamplesFound

        return total_cnt_by_sample

    def general_stats_cols(
        self,
        cnt_by_sample: Dict[str, int],
        pct_by_taxon_by_rank: Dict[str, Dict[str, float]],
        cnt_by_taxon_by_rank_by_sample: Dict[str, Dict[str, Dict[str, int]]],
    ):
        """Add a couple of columns to the General Statistics table"""

        # Get top taxa in most specific taxa rank that we have
        top_taxa = []
        top_rank_code = None
        top_rank_name = None
        for rank_code, rank_name in MultiqcModule.T_RANKS.items():
            if rank_code not in pct_by_taxon_by_rank:
                continue

            sorted_items = sorted(pct_by_taxon_by_rank[rank_code].items(), key=lambda x: x[1], reverse=True)
            for taxon, pct_sum in sorted_items[: MultiqcModule.TOP_N]:
                top_taxa.append(taxon)
            top_rank_code = rank_code
            top_rank_name = rank_name
            break

        if not top_taxa or not top_rank_code or not top_rank_name:
            log.error("No taxa found")
            return

        # Column headers
        headers = dict()

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
        table_pct_by_sample: Dict[str, Dict[str, float]] = {}
        for s_name, cnt_by_taxon_by_rank in cnt_by_taxon_by_rank_by_sample.items():
            _counts = {
                "pct_top_one": cnt_by_taxon_by_rank[top_rank_code].get(top_taxa[0], 0),
                "pct_top_n": sum(cnt_by_taxon_by_rank[top_rank_code].get(t, 0) for t in top_taxa),
            }
            unclassified = cnt_by_taxon_by_rank["U"].get("unclassified")
            if unclassified:
                _counts["pct_unclassified"] = unclassified
            # Convert to percentages
            table_pct_by_sample[s_name] = {k: v / cnt_by_sample[s_name] * 100 for k, v in _counts.items()}

        self.general_stats_addcols(table_pct_by_sample, headers)

    def top_taxa_barplot(
        self,
        total_cnt_by_sample: Dict[str, int],
        pct_by_top_taxon_by_rank: Dict[str, Dict[str, float]],
        cnt_by_top_taxon_by_rank_by_sample: Dict[str, Dict[str, Dict[str, int]]],
    ):
        """Add a bar plot showing the top-N from each taxa rank"""

        rank_datasets = []
        cats = list()
        # Keeping track of encountered codes to display only tabs with available data
        found_rank_codes = set()

        for rank_code in [r for r in MultiqcModule.T_RANKS if r not in ["U", "R"]]:
            if rank_code not in pct_by_top_taxon_by_rank:
                # Taxa rank not found in this dataset
                continue

            rank_cats = dict()
            rank_cnt_data_by_taxon_by_sample: Dict[str, Dict[str, int]] = defaultdict(lambda: defaultdict(int))
            rank_counts_shown: Dict[str, int] = defaultdict(int)
            i = 0
            # Loop through the summed tax percentages to get the top-N across all samples
            for taxon, pct_sum in sorted(pct_by_top_taxon_by_rank[rank_code].items(), key=lambda x: x[1], reverse=True):
                i += 1
                if i > MultiqcModule.TOP_N:
                    break
                rank_cats[taxon] = {"name": taxon}
                # Pull out counts for this rank + classif from each sample
                for s_name, cnt_by_taxon_by_rank in cnt_by_top_taxon_by_rank_by_sample.items():
                    if rank_code not in cnt_by_taxon_by_rank:
                        # Taxa rank not found in this sample
                        continue

                    found_rank_codes.add(rank_code)

                    cnt = cnt_by_taxon_by_rank[rank_code].get(taxon, 0)
                    rank_cnt_data_by_taxon_by_sample[s_name][taxon] += cnt
                    rank_counts_shown[s_name] += cnt

            # Add unclassified count to each rank-level dataset
            for s_name, cnt_by_top_taxon_by_rank in cnt_by_top_taxon_by_rank_by_sample.items():
                cnt = cnt_by_top_taxon_by_rank_by_sample[s_name].get("U", {}).get("unclassified", 0)
                rank_cnt_data_by_taxon_by_sample[s_name]["U"] = cnt
                rank_counts_shown[s_name] += cnt

            # Add in unclassified reads and "other" - we presume from other species etc.
            for s_name, cnt_by_top_taxon_by_rank in cnt_by_top_taxon_by_rank_by_sample.items():
                rank_cnt_data_by_taxon_by_sample[s_name]["other"] = (
                    total_cnt_by_sample[s_name] - rank_counts_shown[s_name]
                )

                # This should never happen... But it does sometimes if the total read count is a bit off
                if rank_cnt_data_by_taxon_by_sample[s_name]["other"] < 0:
                    log.debug(
                        "Found negative 'other' count for {} ({}): {}".format(
                            s_name, MultiqcModule.T_RANKS[rank_code], rank_cnt_data_by_taxon_by_sample[s_name]["other"]
                        )
                    )
                    rank_cnt_data_by_taxon_by_sample[s_name]["other"] = 0

            rank_cats["other"] = {"name": "Other", "color": "#cccccc"}
            rank_cats["U"] = {"name": "Unclassified", "color": "#d4949c"}

            cats.append(rank_cats)
            rank_datasets.append(rank_cnt_data_by_taxon_by_sample)

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
            plot=bargraph.plot(rank_datasets, cats, pconfig),
        )

    def top_taxa_duplication_heatmap(
        self,
        pct_by_top_taxon_by_rank: Dict[str, Dict[str, float]],
        species_minimizer_duplication_by_top_taxon_by_sample: Dict[str, Dict[str, float]],
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

        SPECIES_CODE = "S"
        # Loop through the summed tax percentages to get the top taxa across all samples
        if SPECIES_CODE not in pct_by_top_taxon_by_rank:
            log.debug(f"Taxa rank {SPECIES_CODE} not found, skipping taxa duplication heatmap")
            return

        dup_by_taxon_by_sample: Dict[str, Dict[str, Union[int, None]]] = defaultdict(lambda: defaultdict(int))
        pct_by_top_taxon = pct_by_top_taxon_by_rank[SPECIES_CODE]
        sorted_items = list(sorted(pct_by_top_taxon.items(), key=lambda x: x[1], reverse=True))
        top_sorted_items = sorted_items[: MultiqcModule.TOP_N]
        for taxon, pct_sum in top_sorted_items:
            # Pull out counts for this rank + classif from each sample
            for s_name, dup_by_taxon in species_minimizer_duplication_by_top_taxon_by_sample.items():
                minimizer_duplication = dup_by_taxon.get(taxon, None)
                if minimizer_duplication:
                    dup_by_taxon_by_sample[s_name][taxon] = int(minimizer_duplication)

        # Strip empty samples
        for sample, vals in dict(dup_by_taxon_by_sample).items():
            if len(vals) == 0:
                del dup_by_taxon_by_sample[sample]

        if not dup_by_taxon_by_sample:
            return

        # Build data structures for heatmap
        samples = list(dup_by_taxon_by_sample.keys())
        top_taxa = list(dict(top_sorted_items).keys())
        for sample in dup_by_taxon_by_sample:
            duplication.append(list(dup_by_taxon_by_sample[sample].values()))

        self.add_section(
            name="Duplication rate of top species",
            anchor=f"{self.anchor}-duplication-topfive",
            description=f"The duplication rate of minimizer falling into the top {MultiqcModule.TOP_N} species",
            helptext=f"""
                To make this plot, the minimizer duplication rate is computed for the top {MultiqcModule.TOP_N} most abundant species in all samples.

                The minimizer duplication rate is defined as: `duplication rate = (total number of minimizers / number of distinct minimizers)`

                A low coverage and high duplication rate (`>> 1`) is often sign of read stacking, which probably indicates of false positive hit.
            """,
            plot=heatmap.plot(
                duplication,
                xcats=samples,
                ycats=top_taxa,
                pconfig=pconfig,
            ),
        )


def parse_logs(
    f,
) -> Tuple[
    Dict[str, Dict[str, Union[int]]],
    Dict[str, float],
]:
    """
    Parse a kraken report output file. Only take the top ranks.

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
    (optional, only in new version with minimizers) 7. NCBI taxonomic ID number
    (optional, only in new version with minimizers) 8. Indented scientific name
    """

    cnt_by_rank_by_taxon: Dict[str, Dict[str, int]] = defaultdict(dict)
    min_dup_by_taxon: Dict[str, float] = dict()

    for i, line in enumerate(f["f"]):
        fields = line.split("\t")
        if len(fields) < 6:
            log.error(f"Error parsing Kraken report: {f['fn']} line {i+1} has less than 6 fields: {line}")
            return {}, {}

        if len(fields) == 8:
            # if 8 fields, the new log experimental log (with distinct minimizer)
            (
                percent,
                counts_rooted,
                counts_direct,
                minimizer,
                minimizer_distinct,
                rank_code,
                tax_id,
                taxon,
            ) = fields
        else:
            # If 6 fields are used, it's the 'old' log (without distinct minimizer)
            percent, counts_rooted, counts_direct, rank_code, tax_id, taxon = fields[:6]
            minimizer = None
            minimizer_distinct = None

        taxon = taxon.strip()
        if taxon == "root":
            rank_code = "R"  # can be "-" sometimes

        # This check will skip a lot of lines on the real-life data!
        if rank_code not in MultiqcModule.T_RANKS:
            continue

        counts_rooted = int(counts_rooted)

        # percent = float(percent)
        # counts_direct = int(counts_direct)
        # tax_id = int(tax_id)
        # num_spaces = len(classif) - len(classif_stripped)

        cnt_by_rank_by_taxon[rank_code][taxon] = counts_rooted
        if minimizer_distinct is not None and rank_code == "S":
            minimizer_duplication = int(minimizer) / int(minimizer_distinct) if minimizer_distinct != 0 else 0.0
            min_dup_by_taxon[taxon] = minimizer_duplication

    if "R" not in cnt_by_rank_by_taxon:
        # Can be missing in case if all reads are unassigned
        cnt_by_rank_by_taxon["R"] = {"root": 0}

    return cnt_by_rank_by_taxon, min_dup_by_taxon
