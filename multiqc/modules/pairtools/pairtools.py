import os
from typing import Dict

import yaml
import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph
from multiqc import config

from itertools import zip_longest
from multiqc.modules.pairtools.utils import read_stats_from_file


log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="pairtools",
            anchor="pairtools",
            href="https://github.com/mirnylab/pairtools",
            info="Toolkit for Chromatin Conformation Capture experiments. Handles short-reads paired reference "
            "alignments, extracts 3C-specific information, and perform common tasks such as sorting, filtering, "
            "and deduplication.",
            doi=["10.5281/zenodo.1490830", "10.1101/2023.02.13.528389"],
        )

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        # Find and load any pairtools stats summary files:
        self.pairtools_stats = {}
        for f in self.find_log_files("pairtools", filehandles=True):
            s_name = f["s_name"]
            if (_sample_report := self.parse_pairtools_stats(f)) is None:
                log.warning(f"{s_name} is missing important metrics will not be reported !")
                continue
            else:
                self.pairtools_stats[s_name] = _sample_report
                # Add file to multiqc_sources.txt
                self.add_data_source(f, s_name=s_name)

        # Filter to strip out ignored sample names
        self.pairtools_stats = self.ignore_samples(self.pairtools_stats)
        if len(self.pairtools_stats) == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {len(self.pairtools_stats)} reports")

        # Write parsed report data to a file
        self.write_data_file(self.pairtools_stats, "multiqc_pairtools")

        # Add to self.js to be included in template
        self.js = {
            "assets/js/multiqc_pairtools.js": os.path.join(
                os.path.dirname(__file__), "assets", "js", "multiqc_pairtools.js"
            )
        }

        # load various parameters stored in a separate yml (e.g. color schemes)
        with open(os.path.join(os.path.dirname(__file__), "assets", "params", "params.yml"), "r") as fp:
            self.params = yaml.safe_load(fp)

        self.pairtools_general_stats()

        # Report sections
        self.add_section(
            name="Pairs by alignment status",
            anchor="pair-types",
            description="Number of pairs classified according to their alignment status,"
            " including uniquely mapped (UU), unmapped (NN), duplicated (DD), and others.",
            helptext="""For further details check
            <a href=\"https://pairtools.readthedocs.io/en/latest/formats.html#pair-types\" > pairtools</a>
            documentation.""",
            plot=self.pair_types_chart(),
        )

        self.add_section(
            name="Pre-filtered pairs by genomic location (overview)",
            anchor="cis-ranges-trans",
            description="Distribution of pre-filtered pairs (mapping uniquely and rescued) by genomic"
            " separation for <it>cis-</it>pairs and <it>trans-</it>pairs.",
            helptext="""
            Samples can have different distributions of pairs by genomic locations for various biological
            and technical reasons. Biological examples: cell-cycle stages, differentitation stages difference,
            mutations affecting genome organization, etc.
            Technical differences arise due to the fact that pre-filtered read pairs still include artifacts:
            Short-range cis-pairs (<1kb) are typically enriched in ligation artifacts (self-circles, dangling-ends, etc).
            Elevated number of trans interactions typically suggests increased noise levels - e.g. random ligations etc.""",
            plot=self.pairs_by_cisrange_trans(),
        )

        self.add_section(
            name="Pre-filtered pairs as a function of genomic separation (in detail)",
            anchor="scalings-plots",
            description="Number of interactions (pre-filtered pairs) as a function"
            " of genomic separation, log-binned."
            " Click on an individual curve to reveal information for different"
            " read pair orientations.",
            helptext="""Short-range cis-pairs are typically enriched in ligation artifacts.
            Frequency of interactions for read pairs of different orientations
            ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into these artifacts.
            For example, dangling-ends manifest themselves as FR-pairs, while self-circles - RF.
            Thus, enrichment of FR/RF pairs at a given genomic separation can hint at the level
            of contamination.""",
            plot=self.pairs_with_genomic_separation(),
        )

        self.add_section(
            name="Fraction of read pairs by strand orientation",
            anchor="read-orientation",
            description="Number of pre-filtered pairs reported for every type"
            " of read pair orientation. Numbers are reported for different"
            " ranges of genomic separation and combined.",
            helptext="""Short-range cis-pairs are typically enriched in technical artifacts.
            Frequency of interactions for read pairs of different orientations
            ++,+-,-+ and -- (FF, FR, RF, RR) provide insight into these technical artifacts.
            For example, dangling-ends manifest themselves as FR-pairs, while self-circles - RF.
            Thus, enrichment of FR/RF pairs at a given genomic separation can hint at the level
            of contamination.""",
            plot=self.pairs_by_strand_orientation(),
        )

    @staticmethod
    def parse_pairtools_stats(f):
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
        in a pre-defined order and coloring.
        """

        # Construct a data structure for the plot: { sample: {ptype : count} }
        # and keep track of the observed pair types - for consistent visuals:
        data_dict = {}
        _observed_ptypes = set()
        for _sample_name, _sample_data in self.pairtools_stats.items():
            if (_sample_data_dict := _sample_data["pair_types"]) is None:
                return None  # if any of samples is None -> None
            else:
                data_dict[_sample_name] = _sample_data_dict
                _observed_ptypes |= set(_sample_data_dict)

        # display common pairtypes first with pre-defined colors and order :
        ptypes_annotated: Dict = {}
        for _ptype, ptype_color in self.params["pairtypes_colors"].items():
            if any(t.lower() == _ptype.lower() for t in _observed_ptypes):
                ptypes_annotated[_ptype] = {"color": ptype_color, "name": _ptype}
            else:
                # display remaining pairtypes second with automatic colors :
                ptypes_annotated[_ptype] = {"name": _ptype}

        # multiqc interactive plotting call :
        return bargraph.plot(
            data_dict,
            ptypes_annotated,
            pconfig=bargraph.BarPlotConfig(
                id="pair_types",
                title="pairtools: pair types report",
                ylab="# Reads",
                cpswitch_counts_label="Number of Reads",
            ),
        )

    def pairs_by_cisrange_trans(self):
        """
        cis-pairs split into several ranges
        of genomic separation, and trans-category.

        Here we extract relevant sub-dict "cis_dist" for
        each sample, make sure there are no None-s and
        dist_names are the same across all samples.
        """

        data_dict = {}
        _dist_category_set = set()
        for _sample_name, _sample_data in self.pairtools_stats.items():
            if (_sample_data_dict := _sample_data["cis_dist"]) is None:
                return None  # if any of samples is None -> None
            else:
                data_dict[_sample_name] = _sample_data_dict
                # make sure all dist_categories are identical
                _dist_category_set |= _sample_data_dict.keys()
                if len(_dist_category_set) > 1:
                    return None

        # extract uniq distance categories
        uniq_dist_categories = _dist_category_set

        # assign colors to distance ranges (as many as possible):
        key_dict: Dict = {}
        for color, _cat_names in zip_longest(self.params["cis_range_colors"], uniq_dist_categories):
            for _cat_name in _cat_names:
                key_dict[_cat_name] = {"name": _cat_name}
                if color:
                    key_dict[_cat_name]["color"] = color
                if _cat_name == "trans":
                    key_dict[_cat_name]["color"] = "#1035AC"  # dark-ish blue

        # multiqc interactive plotting call :
        return bargraph.plot(
            data_dict,
            key_dict,
            pconfig=bargraph.BarPlotConfig(
                id="pair_cis_ranges",
                title="pairtools: cis pairs broken into ranges",
                ylab="# Reads",
                cpswitch_counts_label="Number of reads",
            ),
        )

    def pairs_by_strand_orientation(self):
        """
        number of cis-pairs with genomic separation
        """

        # traverse sample_data to extract dist_categories
        _dist_categories = None
        for _sample_name, _sample_data in self.pairtools_stats.items():
            if (_sample_data_dict := _sample_data["pairs_by_strand"]) is None:
                return None  # if any of samples is None -> None
            else:
                # extract distance categories that should be identical across sample
                _dist_categories = list(_sample_data_dict.keys())

        # initialize data_dict - with distance categories as the outter layer
        if _dist_categories is None:
            return None

        data_dict: Dict[str, Dict] = {_cat: {} for _cat in _dist_categories}

        # traverse sample_data again to reorder data for plotting
        for _sample_name, _sample_data in self.pairtools_stats.items():
            _sample_data_dict = _sample_data["pairs_by_strand"]
            for _cat, _innermost_dict in _sample_data_dict.items():
                data_dict[_cat][_sample_name] = _innermost_dict

        # annotate read orientations with nice colors:
        keys_annotated = {}
        for key in ["++", "-+", "+-", "--"]:
            name = self.params["pairs_orientation_names"][key]
            color = self.params["pairs_orientation_colors"][key]
            keys_annotated[key] = {"color": color, "name": name}

        return bargraph.plot(
            [data_dict[_dist] for _dist in _dist_categories],
            [keys_annotated for _dist in _dist_categories],
            pconfig=bargraph.BarPlotConfig(
                id="pair_by_orient_cis_ranges",
                title="Cis pairs broken into ranges and read orientations",
                ylab="# Reads",
                cpswitch_counts_label="Number of reads",
                data_labels=list(data_dict.keys()),
            ),
        )

    def pairs_with_genomic_separation(self):
        """
        number of cis-pairs with genomic separation
        # dist_freq/56234133-100000000/-+
        """

        _strand_categories = ["all", "++", "+-", "--", "-+"]
        # initialize data_dict - with distance categories as the outter layer
        data_dict: Dict[str, Dict] = {_cat: {} for _cat in _strand_categories}

        # traverse sample_data to extract dist_categories
        for _sample_name, _sample_data in self.pairtools_stats.items():
            if (_sample_data_dict := _sample_data["dist_freq"]) is None:
                return None  # if any of samples is None -> None
            # reorder data in _sample_data_dict to extract strand_category on top:
            for _cat in _strand_categories:
                data_dict[_cat][_sample_name] = _sample_data_dict[_cat]

        return linegraph.plot(
            [data_dict[_cat] for _cat in _strand_categories],
            pconfig=linegraph.LinePlotConfig(
                id="broom_plot",
                title="pairtools: Pairs by distance and orientation",
                xlab="Genomic separation (bp)",
                ylab="number of pairs",
                xlog=True,
                ylog=True,
                data_labels=[
                    {"name": "P(s)", "ylab": "number of pairs"},
                    {"name": "FF", "ylab": "number of pairs"},
                    {"name": "FR", "ylab": "number of pairs"},
                    {"name": "RF", "ylab": "number of pairs"},
                    {"name": "RR", "ylab": "number of pairs"},
                ],
            ),
        )

    def pairtools_general_stats(self):
        """Add columns to General Statistics table"""
        headers = {
            "total": {
                "title": "Read pairs",
                "description": "Total read pairs before mapping",
                "scale": "Greys",
                "shared_key": "read_count",
            },
            "frac_unmapped": {
                "title": "Unmapped",
                "description": "% of pairs (w.r.t. total) with both sides unmapped",
                "suffix": "%",
                "scale": "PuRd",
            },
            "frac_single_sided_mapped": {
                "title": "One-sided",
                "description": "% of pairs (w.r.t. total) with one side mapped",
                "suffix": "%",
                "scale": "Purples",
            },
            "frac_mapped": {
                "title": "Two-sided",
                "description": "% of pairs (w.r.t. total) with both sides mapped",
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "frac_dups": {
                "title": "Duplicated",
                "description": "% of duplicated pairs (w.r.t. mapped)",
                "suffix": "%",
                "scale": "YlOrRd",
            },
            "total_nodups": {
                "title": "Unique pairs",
                "description": "Mapped pairs after deduplication",
                "scale": "Greys",
                "shared_key": "read_count",
            },
            "frac_cis": {
                "title": "Cis",
                "description": "% of cis-pairs (w.r.t mapped)",
                "suffix": "%",
                "scale": "PuOr-rev",
            },
        }

        self.general_stats_addcols(self.pairtools_stats, headers)