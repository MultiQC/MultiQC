#!/usr/bin/env python

""" MultiQC module to parse stats output from pairtools """

import logging
import os
from collections import OrderedDict
from copy import copy
from itertools import zip_longest
from random import choice

import numpy as np
import yaml

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, heatmap, linegraph

from .utils import (
    contact_areas_genomewide,
    cumsums_to_rangesums,
    edges_to_intervals,
    genomic_dist_human_str,
    read_stats_from_file,
    total_coverage,
)

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    """This MultiQC module parses various
    stats produced by pairtools."""

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="pairtools",
            anchor="pairtools",
            href="https://github.com/mirnylab/pairtools",
            info="pairtools is a command-line framework for processing sequencing data"
            " generated with Chromatin Conformation Capture based experiments:"
            " pairtools can handle pairs of short-reads aligned to a reference genome,"
            " extract 3C-specific information and perform common tasks, such as sorting,"
            " filtering and deduplication.",
            doi="10.5281/zenodo.1490831",
        )

        # Find and load any pairtools stats summary files:
        self.pairtools_stats = dict()
        for f in self.find_log_files("pairtools", filehandles=True):
            s_name = f["s_name"]
            self.pairtools_stats[s_name] = self.parse_pairtools_stats(f)

        # Filter to strip out ignored sample names
        self.pairtools_stats = self.ignore_samples(self.pairtools_stats)

        if len(self.pairtools_stats) == 0:
            raise UserWarning("No reports to use.")

        log.info(f"Found {len(self.pairtools_stats)} reports")

        # Add to self.js to be included in template
        self.js = {
            "assets/js/multiqc_pairtools.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "multiqc_pairtools.js"
            )
        }

        # load various parameters stored in a separate yml (e.g. color schemes)
        with open(os.path.join(os.path.dirname(__file__), "assets", "params", "params.yml"), "r") as fp:
            self.params = yaml.safe_load(fp)

        #############################################
        #############################################
        # should add several cross-sample checks
        #############################################
        #############################################
        # e.g. check if all sets of distances are identical,
        # by comparing them all to the last one: (s_name, sorted_dists)
        chromsizes_provided = True
        if len(self.pairtools_stats) > 1:
            random_sample = choice(list(self.pairtools_stats))
            # check if cis-ranges are same across samples
            random_dists = self.pairtools_stats[random_sample]["cis_dist"]["dists"]
            for s_name in self.pairtools_stats:
                # check is cis_dists (e.g. cis_10kb+) have identical dists across samples
                # it is [1, 2, 4, 10, 20, 40] as of now in pairtools
                if self.pairtools_stats[s_name]["cis_dist"]["dists"] != random_dists:
                    log.warning(
                        f"Samples {s_name} and {random_sample} have different sets of cis-ranges,\n"
                        "pairs by cis range will not be reported !"
                    )
            # check if chromsizes are the same across samples
            random_chromsizes = self.pairtools_stats[random_sample]["chromsizes"]
            for s_name in self.pairtools_stats:
                if not self.pairtools_stats[s_name]["chromsizes"]:
                    chromsizes_provided = False
                elif self.pairtools_stats[s_name]["chromsizes"] != random_chromsizes:
                    chromsizes_provided = False
                    log.warning(f"Samples {s_name} and {random_sample} have different sets of chromsizes.")

        # determine max total reads for general stats:
        self.max_total_reads = 0
        for s_name in self.pairtools_stats:
            self.max_total_reads = max(self.max_total_reads, self.pairtools_stats[s_name]["total"])

        self.pairtools_general_stats()

        # Report sections
        self.add_section(
            name="Pairs alignment status",
            anchor="pair-types",
            description="Number of pairs classified according to their alignment status,"
            " including uniquely mapped (UU), unmapped (NN), duplicated (DD), and others.",
            helptext="""For further details check
                        <a href=\"https://pairtools.readthedocs.io/en/latest/formats.html#pair-types\" > pairtools</a>
                        documentation.""",
            plot=self.pair_types_chart(),
        )

        self.add_section(
            name="Pre-filtered pairs grouped by genomic separations",
            anchor="cis-ranges-trans",
            description="Distribution of pre-filtered pairs (UU, UR and RU) by genomic"
            " separation for <it>cis-</it>pairs and <it>trans-</it>pairs.",
            helptext="""Pre-filtered read pairs might still include artifacts:
            Short-range cis-pairs (<1kb) are typically enriched in technical artifacts (self-circles, dangling-ends, etc).
            High fraction of trans interactions typically suggests increased noise levels""",
            plot=self.pairs_by_cisrange_trans(),
        )

        if not chromsizes_provided:
            nochromsizes_warning = """**Warning! There are no chromosome sizes detected in the stats files, therefore
            interaction frequencies are presented without normalization (i.e. raw log-binned counts)**"""
        else:
            nochromsizes_warning = ""
        self.add_section(
            name="Frequency of interactions as a function of genomic separation",
            anchor="scalings-plots",
            description="Frequency of interactions (pre-filtered pairs) as a function"
            ' of genomic separation, known as "scaling plots", P(s).'
            " Click on an individual curve to reveal P(s) for different"
            " read pair orientations." + nochromsizes_warning,
            helptext="""Short-range cis-pairs are typically enriched in technical artifacts.
            Frequency of interactions for read pairs of different orientations
            ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into these technical artifacts.
            For example, dangling-ends manifest themselves as FR-pairs, while self-circles - RF.
            Thus enrichment of FR/RF pairs at a given genomic separation can hint at the level
            of contamination.""",
            plot=self.pairs_with_genomic_separation(),
        )

        self.add_section(
            name="Fraction of read pairs by strand orientation",
            anchor="read-orientation",
            description="Number of interactions (pre-filtered pairs) reported for every type"
            " of read pair orientation. Numbers are reported for different"
            " ranges of genomic separation and combined.",
            helptext="""Short-range cis-pairs are typically enriched in technical artifacts.
            Frequency of interactions for read pairs of different orientations
            ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into these technical artifacts.
            For example, dangling-ends manifest themselves as FR-pairs, while self-circles - RF.
            Thus enrichment of FR/RF pairs at a given genomic separation can hint at the level
            of contamination.""",
            plot=self.pairs_by_strand_orientation(),
        )

        if chromsizes_provided:
            self.add_section(
                name="Pre-filtered pairs grouped by chromosomes",
                anchor="pairs-by-chroms",
                description="Number of pre-filtered interactions (pairs) within a single chromosome"
                " or for a pair of chromosomes.",
                helptext="""Numbers of pairs are normalized by the total number of pre-filtered pairs per sample.
                Number are reported only for chromosomes/pairs that have >1% of pre-filtered pairs.""",
                plot=self.coverage_by_chrom(),
            )

    def parse_pairtools_stats(self, f):
        """
        Parse a pairtools summary stats
        """
        f_handle = f["f"]
        log.info(f"parsing .stats file: {f_handle.name}")
        return read_stats_from_file(f_handle)

    def pair_types_chart(self):
        """
        Generate the pair types report: a stacked bargraph
        with the number of reads per pair type displayed
        in a pre-defined orders and a set of colors.
        """
        ptypes_field = "pair_types"

        # Construct a data structure for the plot: { sample: {ptype : count} }
        # and keep track of the observed pair types - for nice visuals:
        ptypes_dict = dict()
        observed_ptypes = set()
        for s_name in self.pairtools_stats:
            ptypes_dict[s_name] = copy(self.pairtools_stats[s_name][ptypes_field])
            observed_ptypes |= set(ptypes_dict[s_name])

        # display common pairtypes first with pre-defined colors and order :
        ptypes_anotated = OrderedDict()
        for common_ptype, ptype_color in self.params["pairtypes_colors"].items():
            if common_ptype in observed_ptypes:
                ptypes_anotated[common_ptype] = {"color": ptype_color, "name": common_ptype}
                observed_ptypes.discard(common_ptype)
        # display remaining pairtypes second with automatic colors :
        for rare_ptype in observed_ptypes:
            ptypes_anotated[common_ptype] = {"name": rare_ptype}

        # config for the pairtype barplot :
        config = {
            "id": ptypes_field,
            "title": "pairtools: pair types report",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        # multiqc interactive plotting call :
        return bargraph.plot(ptypes_dict, ptypes_anotated, pconfig=config)

    def pairs_by_cisrange_trans(self):
        """
        cis-pairs split into several ranges
        of genomic separation, and trans-category.

        Several important assumptions made here:
         - provided counts are cumulative
         - "distance" intervals are "oppen", e.g. (1kb+,5kb+,...)
         - genomic "distances" are in kb
        """

        # Construct a data structure for the plot
        cis_rangecounts_dict = dict()
        cis_distances_dict = dict()
        for s_name in self.pairtools_stats:
            sample_stats = self.pairtools_stats[s_name]
            sorted_dists = sample_stats["cis_dist"]["dists"]  # edges of distance ranges
            cumcounts = sample_stats["cis_dist"]["counts"]  # cumulative counts
            rangecounts = cumsums_to_rangesums(cumcounts, sample_stats["cis"])
            rangecounts.append(sample_stats["trans"])
            cis_rangecounts_dict[s_name] = rangecounts
            cis_distances_dict[s_name] = sorted_dists

        # check if all sets of distances are identical,
        # by comparing them all to the last one: (s_name, sorted_dists)
        for _s, _d in cis_distances_dict.items():
            if _d != sorted_dists:
                log.warning(
                    f"Samples {s_name} and {_s} have different sets of cis-ranges,\n"
                    "pairs by cis range will not be reported !"
                )
                # it is [1, 2, 4, 10, 20, 40] as of now in pairtools, but will it stay like that ?
                # return bargraph.plot({}, {}, )  # returning empty barplot
                return None

        # generate nice distance ranges for printing :
        range_names = []
        for start, end in edges_to_intervals(sorted_dists):
            if start == 0:  # shortest cis range
                range_names.append(f"cis: <{end}kb")
            elif end is None:  # longest cis range
                range_names.append(f"cis: >{start}kb")
            else:  # all in betweens
                range_names.append(f"cis: {start}-{end}kb")
        # trans pairs (ultimate long range):
        range_names.append("trans")
        # now zip counts and range names together for plotting:
        data_dict = {s: dict(zip(range_names, c)) for s, c in cis_rangecounts_dict.items()}

        # color distance ranges with nice pre-defined colors (as many as possible):
        key_dict = OrderedDict()
        for color, dist_range in zip_longest(self.params["cis_range_colors"], range_names):
            key_dict[dist_range] = {"name": dist_range}
            if color:
                key_dict[dist_range]["color"] = color

        # Config for the plot
        config = {
            "id": "pair_cis_ranges",
            "title": "pairtools: cis pairs broken into ranges",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        # multiqc interactive plotting call :
        return bargraph.plot(data_dict, key_dict, pconfig=config)

    def pairs_by_strand_orientation(self):
        """
        number of cis-pairs with genomic separation
        """

        dist_freq_field = "dist_freq"
        orientation_keys = ["++", "-+", "+-", "--"]
        dist_ranges = edges_to_intervals(self.params["pairs_distance_edges"])

        # Construct a data structure for the plot
        _data = dict()
        for key in dist_ranges:
            _data[key] = dict()
            for s_name in self.pairtools_stats:
                _data[key][s_name] = dict()

        for s_name in self.pairtools_stats:
            # extract scalings data structure per sample:
            sample_bins = self.pairtools_stats[s_name]["dist_bins"]
            sample_dist_freq = self.pairtools_stats[s_name][dist_freq_field]
            # given a set of fixed distance ranges extract FF,RR,FR,RF:
            for start, end in dist_ranges:
                # determine start index, according to "dist_bins"
                start_idx = np.searchsorted(sample_bins, start) if start else 0
                # determine end index, according to "dist_bins"
                end_idx = np.searchsorted(sample_bins, end) if end is not None else None
                # calculate ratios of FF FR RF RR ...
                for orient in orientation_keys:
                    # slice of the scaling for a given range of distances:
                    sliced_data = sample_dist_freq[orient][start_idx:end_idx]
                    _data[(start, end)][s_name][orient] = np.sum(sliced_data.astype(float))

        # generate pretty labels for ranges of genomic separation
        data_labels = []
        for start, end in dist_ranges:
            start_str = genomic_dist_human_str(start)
            end_str = genomic_dist_human_str(end)
            if start == 0:
                data_labels.append(f"<{end_str}")
            elif end is None:
                data_labels.append(f">{start_str}")
            else:
                data_labels.append(f"{start_str}-{end_str}")

        # Config for the plot
        config = {
            "id": "pair_by_orient_cis_ranges",
            "title": "pairtools: cis pairs broken into ranges and read orientations",
            "ylab": "# Reads",
            "cpswitch_counts_label": "Number of Reads",
            "data_labels": data_labels,
        }

        # annotate read orientations with nice colors:
        keys_annotated = OrderedDict()
        for key in orientation_keys:
            name = self.params["pairs_orientation_names"][key]
            color = self.params["pairs_orientation_colors"][key]
            keys_annotated[key] = {"color": color, "name": name}

        return bargraph.plot(
            [_data[_d] for _d in dist_ranges],
            [keys_annotated for _d in dist_ranges],
            pconfig=config,
        )

    def pairs_with_genomic_separation(self):
        """
        number of cis-pairs with genomic separation
        # dist_freq/56234133-100000000/-+
        """

        dist_freq_field = "dist_freq"
        data_cats = ["avg", "FF", "FR", "RR", "RF"]
        # extract pair orientations types - as they'll be used here many times
        porient_names = self.params["pairs_orientation_names"]

        # Initialize a data structure for the plot
        _data = {cat: {} for cat in data_cats}

        for s_name in self.pairtools_stats:
            for cat in data_cats:
                _data[cat][s_name] = {}
            # pre-calculate geom-mean of dist-bins for P(s):
            _dist_bins = np.asarray(self.pairtools_stats[s_name]["dist_bins"])
            # extract contacts by distance ranges and orientation types
            sample_dist_freq = self.pairtools_stats[s_name][dist_freq_field]

            # contacts by orientation and distance, and their average summary,
            _summary = {}
            _summary["avg"] = np.sum([sample_dist_freq[po] for po in porient_names], axis=0).astype(float)
            for po, po_name in porient_names.items():
                _summary[po_name] = np.asarray(sample_dist_freq[po]).astype(float)

            if self.pairtools_stats[s_name]["chromsizes"]:
                # normalize contacts by distance with the theoretical # of pairs in a range
                sizes = np.fromiter(self.pairtools_stats[s_name]["chromsizes"].values(), dtype=float)
                _areas = contact_areas_genomewide(_dist_bins, scaffold_sizes=sizes)
                for cat in data_cats:
                    _summary[cat] = _summary[cat] / _areas

            # assign geometric mean distance to every distance interval
            _dist_bins_geom = np.sqrt(_dist_bins[:-1] * _dist_bins[1:])
            # fill in the data for XY-line plotting
            # i.e. dict by samples of dicts by dist (X), of (normalized) counts (Y):
            for cat in data_cats:
                _data[cat][s_name] = dict(zip(_dist_bins_geom[1:-1], _summary[cat][1:]))

        pconfig = {
            "id": "broom_plot",
            "title": "Pairs by distance and by read orientation",
            "xlab": "Genomic separation (bp)",
            "xLog": True,
            "yLog": True,
            "data_labels": [
                {"name": "P(s)", "ylab": "frequency of interactions"},
                {"name": "FF", "ylab": "frequency of interactions"},
                {"name": "FR", "ylab": "frequency of interactions"},
                {"name": "RF", "ylab": "frequency of interactions"},
                {"name": "RR", "ylab": "frequency of interactions"},
            ],
            "click_func": "single_scaling",  # activate custom JS when individual curve clicked
        }

        return linegraph.plot([_data[cat] for cat in data_cats], pconfig=pconfig)

    # coverage per chromosome per sample - a heatmap
    def coverage_by_chrom(self):
        """
        number of pairs by chromosome pairs
        """

        the_data = []
        for s_name in self.pairtools_stats:
            # output not more than 100 chroms, sorted on size ...
            if len(self.pairtools_stats[s_name]["chromsizes"]):
                all_chroms = list(self.pairtools_stats[s_name]["chromsizes"])
                sizes = list(self.pairtools_stats[s_name]["chromsizes"].values())
            else:
                # returning empty
                return heatmap.plot([])
            # show first 100 chroms only ...
            cov_chroms = all_chroms[:100]

            tot_contact = self.pairtools_stats[s_name]["total_nodups"] / self.max_total_reads
            coverage = total_coverage(self.pairtools_stats[s_name]["chrom_freq"], cov_chroms)
            # normalize coverage by chromsizes and total nodup interaction per sample
            coverage = [c / (s * tot_contact) for c, s in zip(coverage, sizes)]

            # [ [ _data[s][k] for k in xcats ] for s in ycats ]
            the_data.append(coverage)

        pconfig = {
            "square": False,
            "xcats_samples": False,
            "ycats_samples": True,
        }

        return heatmap.plot(
            the_data,
            xcats=cov_chroms,
            ycats=list(self.pairtools_stats),
            pconfig=pconfig,
        )

    def pairtools_general_stats(self):
        """Add columns to General Statistics table"""
        headers = OrderedDict()
        headers["total"] = {
            "title": f"{config.read_count_prefix} read pairs",
            "description": f"Total read pairs ({config.read_count_desc})",
            "min": 0,
            "modify": lambda x: x * config.read_count_multiplier,
            "scale": "Blues",
        }
        headers["frac_unmapped"] = {
            "title": "% unmapped",
            "description": "% of pairs (w.r.t. total) with both sides unmapped",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "OrRd",
        }
        headers["frac_single_sided_mapped"] = {
            "title": "% single-side mapped",
            "description": "% of pairs (w.r.t. total) with one side mapped",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
        }
        headers["frac_mapped"] = {
            "title": "% both-side mapped",
            "description": "% of pairs (w.r.t. total) with both sides mapped",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "RdYlGn",
        }
        headers["frac_dups"] = {
            "title": "% duplicated",
            "description": "% of duplicated pairs (w.r.t. mapped)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "OrRd",
        }
        headers["cis_percent"] = {
            "title": "% cis",
            "description": "% of cis-pairs (w.r.t mapped)",
            "max": 100,
            "min": 0,
            "suffix": "%",
            "scale": "YlGn",
        }
        self.general_stats_addcols(self.pairtools_stats, headers, "pairtools")
