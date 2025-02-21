import itertools
import json
import logging
import os
import re
from collections import defaultdict
from typing import Dict

import numpy as np

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import bargraph, linegraph, table

log = logging.getLogger(__name__)

# A 64 color palette defined in four gradient blocks. Red->Yellow, Yellow->Green, Green->Blue, Blue->Deep Blue
# The method color_picker uses these to sample a gradient with varying sparseness

# red to yellow
ry_1 = [
    "#EF476F",
    "#F0506E",
    "#F1586E",
    "#F2616D",
    "#F36A6D",
    "#F4726C",
    "#F57B6C",
    "#F6836B",
    "#F78C6B",
    "#F8956A",
    "#F99D69",
    "#FAA669",
    "#FBAF68",
    "#FCB768",
    "#FDC067",
    "#FEC867",
    "#FFD166",
]

# yellow to green
yg_2 = [
    "#FFD166",
    "#EFD16A",
    "#E0D26D",
    "#D0D271",
    "#C1D275",
    "#B1D378",
    "#A2D37C",
    "#92D37F",
    "#83D483",
    "#73D487",
    "#63D48A",
    "#54D48E",
    "#44D592",
    "#35D595",
    "#25D599",
    "#16D69C",
    "#06D6A0",
]

# green to blue
gb_3 = [
    "#06D6A0",
    "#07D1A1",
    "#07CDA2",
    "#08C8A3",
    "#09C3A5",
    "#09BEA6",
    "#0ABAA7",
    "#0BB5A8",
    "#0CB0A9",
    "#0CABAA",
    "#0DA7AB",
    "#0EA2AC",
    "#0E9DAE",
    "#0F98AF",
    "#1094B0",
    "#108FB1",
    "#118AB2",
]

# blue to dark blue
bd_4 = [
    "#118AB2",
    "#1085AB",
    "#107FA4",
    "#0F7A9E",
    "#0E7597",
    "#0E7090",
    "#0D6A89",
    "#0C6582",
    "#0C607C",
    "#0B5B75",
    "#0A556E",
    "#0A5067",
    "#094B60",
    "#08465A",
    "#084053",
    "#073B4C",
]

# the 64-color gradient
grad64 = ry_1[:-1] + yg_2[:-1] + gb_3[:-1] + bd_4

# 9-color transition related to base colors above.
pal_8 = ["#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3"]
rev_8 = pal_8[::-1]


def color_picker(degen):
    """
    Select colours from a 64 colour gradient, where an effort is made to use the full
    width of the gradient depending on the number of elements required.
    :param degen: a list of degeneracies (number of colours per object)
    :return: a flat list of colours equal the sum of degen
    """
    if len(degen) == 1:
        # a single non-ambiguous enzyme, lets make this blue
        if degen[0] not in {1, 4, 16}:
            raise ValueError(f"got {degen[0]} junc_degen values can only be 1, 4 or 16")
        return grad64[0 :: 64 // degen[0]]
    else:
        cols = []
        for n, jd in enumerate(degen):
            if jd not in {1, 4, 16}:
                raise ValueError(f"got {jd} when junc_degen values can only be 1, 4 or 16")
            cols += grad64[16 * n : 16 * n + 16 : 16 // jd]
        return cols


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="qc3C",
            anchor="qc3C",
            href="http://github.com/cerebis/qc3C",
            info="Reference-free and BAM based quality control for Hi-C data",
            extra="""
            qc3C allows researchers to assess the fraction of read-pairs within a Hi-C library that are a product 
            of proximity ligation -- in effect the Hi-C signal strength. Using a k-mer based approach, signal strength 
            is inferred directly from reads and therefore no reference is required. Reference based assessment is also 
            available and can provide further details.
        
            With this information in hand, researchers are able to decide how much sequencing will be needed to achieve 
            their experimental aims.
            """,
            doi="10.1371/journal.pcbi.1008839",
        )

        self.qc3c_data: Dict[str, Dict] = defaultdict(dict)
        # additional members for conditional plotting of per-genotype junction frequency
        self.do_digest_plot = False
        self.digest_junctions: Dict[str, Dict] = defaultdict(dict)

        for f in self.find_log_files("qc3C", filehandles=True):
            self.parse_qc3c_log(f)

        for k in self.qc3c_data:
            self.qc3c_data[k] = self.ignore_samples(self.qc3c_data[k])

        # check that any non-empty records were stored under one of the two analysis modes
        n_reports = len(self.qc3c_data["kmer"]) + len(self.qc3c_data["bam"])
        if n_reports == 0:
            raise ModuleNoSamplesFound

        log.info(f"Found {n_reports} reports")

        if len(self.qc3c_data["bam"]) > 0:
            self.write_data_file(self.qc3c_data["bam"], "multiqc_qc3c_bam")
            log.debug(f"Found {len(self.qc3c_data['bam'])} BAM analysis reports")

            self.add_section(
                name="BAM mode analysis details",
                anchor="qc3C-bam-runtime-parameters",
                description="""
                    This table details various alignment features which are potentially of interest to
                    researchers attempting to assess the quality of a Hi-C library.
                """,
                plot=self.bam_runtime_table(),
            )

            # self.add_section(
            #     name = 'BAM mode read-pair breakdown',
            #     anchor = 'qc3C-bam-hic-fraction',
            #     description = """
            #         This table details various alignment features which are potentially of interest to
            #         researchers attempting to assess the quality of a Hi-C library.
            #     """,
            #     helptext = """
            #         Here, the **adjusted fraction of read-through** events can be taken as a estimate of the fraction of
            #         Hi-C read-pairs within the library.
            #
            #         * _trans_ (inter-reference) and _cis_ (intra-reference) pairs.
            #         * The fraction of _cis_ pairs whose separation is less than 1000 bp.
            #         * The fraction of total insert extent which was unobserable.
            #         * The fraction of accepted reads whose alignment began with a cutsite.
            #         * The fraction of accepted reads whose alignment ended in a cutsite.
            #         * The fraction of accepted reads which fully aligned and ended in a cutsite.
            #         * The fraction of reads which contained a suspected read-through event (across a junction)
            #         * The fraction of read-through events whose 3-prime end further aligned elsewhere.
            #         * The fraction of read-through events adjusted for unobservable extent.
            #     """,
            #     plot = self.bam_signal_table()
            # )

            self.add_section(
                name="BAM mode read parsing",
                anchor="qc3C-bam-acceptance-plot",
                description="""
                    This figure displays a breakdown of proportion of parsed reads rejected due to various
                    criteria and the proportion that were accepted.
                """,
                plot=self.bam_acceptance_plot(),
            )

            # the following table is disabled to simplify report, as it is redundant consider the plot that follows.
            #
            # self.add_section(
            #     name = 'BAM mode HiCPro categories',
            #     anchor = 'qc3C-bam-hicpro-table',
            #     description = """
            #         This table details the read-pair categorisation devised by
            #         [HiC-Pro](https://github.com/nservant/HiC-Pro).
            #     """,
            #     helptext = """
            #         As the field has moved from 6-cutter to 4-cutter enzymes, and subsequently dual-enzyme digests,
            #         the higher density of sites has made this framework less useful, since it has become increasingly
            #         easy to satisfy the intervening site criteria.
            #     """,
            #     plot = self.bam_hicpro_table()
            # )

            self.add_section(
                name="BAM mode HiC-Pro validation",
                anchor="qc3C-bam-valid-plot",
                description="""
                    A visualisation of the read-pair categories devised by
                    [HiC-Pro](https://github.com/nservant/HiC-Pro).
                """,
                helptext="""
                    As the field has moved from 6-cutter to 4-cutter enzymes, and subsequently dual-enzyme digests, the
                    higher density of sites has made this framework less useful, since it has become increasingly easy
                    to satisfy the intervening site criteria.
                """,
                plot=self.bam_valid_plot(),
            )

            self.add_section(
                name="BAM mode long-range pairs",
                anchor="qc3C-bam-longrange-plot",
                description="This plot visualises the breakdown of read-pairs based on separation distance.",
                helptext="""
                    The breakdown of separation distance is only calculated for *cis*-mapping pairs.

                    Ideally, Hi-C proximity ligation should produce many pairs which are greater than 1000 bp apart.
                    However, these statistics are strongly influenced by the state of the reference. For draft assemblies
                    the distance at which pairs can map is limited by the degree of fragmentation and length of contigs.
                    As a result, many more pairs will be categorised as *trans*-mapping and pairs which are truly
                    inter-molecular cannot be distinguished from those which are merely inter-contig.
                """,
                plot=self.bam_longrange_plot(),
            )

            self.add_section(
                name="BAM mode distribution of fragment separation",
                anchor="qc3C-bam-fragment-histogram",
                description="""
                    This figure displays the a normalised histogram of read-pair separation binned uniformly in
                    log-space.

                    Due to the binning strategy, the x-axis is log-scaled and visually accommodates pair separations up
                    to 1 million bp. The inferred insert size for each library is represented by a dashed, grey vertical
                    line. The y-axis is log-scaled by default, allowing the density attributed to long-range pairs to be
                    more easily seen.
                """,
                helptext="""
                    A characteristic of Hi-C libraries, is the presence of a large peak below 1000 bp. qc3C attributes
                    this to regular (and undesirable) shotgun pairs creeping through the Hi-C protocol. The peak is used
                    by qc3C to infer the insert size, which is later employed to estimate unobservable extent of
                    inserts.

                    **Note:** the inferred insert size can be significantly smaller than what a sequencing facility
                    might report the experimentally determined insert size to be. This discrepancy can be explained by
                    the failure to account for the additional adapter sequence when fragments are assessed during
                    library preparation.
                """,
                plot=self.bam_fragment_histogram(),
            )

            if self.do_digest_plot:
                self.add_section(
                    name="BAM mode junction breakdown",
                    anchor="qc3C-bam-junction-plot",
                    description="""
                        This figure displays the frequency at which a library's possible junction sequences
                        are actually observed in the reads. (_Trivial single-digests are ignored_)
                    """,
                    helptext="""
                        For trivial single-enzyme digests, there is only one possible junction sequence and so the
                        result for these experiments are not plotted. For dual-enzyme (such as Phase Genomics) there are
                        four potential junctions, while for dual-enzyme digests with one ambiguous site (such as Arima
                        Genomics) there are 16 possible junction sequences.

                        How efficiently the more complicated library protocols are at producing hybrid junctions is
                        possibly just a point of interest.

                        Junctions are named for which enzymes was responsible for creating the 5' and 3' ends.
                        E.g. `Sau3AI/MluCI` would involve two different enzymes, while `Sau3AI/Sau3AI` only one, as
                        would be the case in a single-enzyme digest. Proceeding the name is the actual junction
                        sequence.

                        The junctions are grouped by their 5' and then 3' enzyme, while the color spectrum used across
                        each bar aims to emphasise these enzymatic sources.

                        **Note:** in BAM mode, the counts **are** controlled for false positives, in the sense that read
                        alignments must terminate at a cutsite, but the read sequence must continue and contain the
                        observed junction.
                    """,
                    plot=self.bam_junction_plot(),
                )

        if len(self.qc3c_data["kmer"]) > 0:
            self.write_data_file(self.qc3c_data["kmer"], "multiqc_qc3c_kmer")
            log.debug(f"Found {len(self.qc3c_data['kmer'])} k-mer analysis reports")

            self.add_section(
                name="K-mer mode runtime details",
                anchor="qc3C-kmer-runtime-parameters",
                description="""
                    This table includes user specified input options, observed read-length and unobservable fraction.
                """,
                plot=self.kmer_runtime_table(),
            )

            self.add_section(
                name="K-mer mode Hi-C fraction",
                anchor="qc3C-kmer-hic-fraction",
                description="This table lists the inferred proportion of Hi-C proximity ligation fragments.",
                helptext="""
                    Here, **Mean adjusted Hi-C fraction** represents the best estimate of the proportion of a library's
                    read-pairs which are a product of proximity ligation. This figure is arrived at by correcting the
                    raw estimate for the fraction of insert extent which was not observable.

                    The observable extent is limited by the length of reads relative to the supplied insert size, as
                    well as a further constraint on flanking sequence around any suspected junction sequence.
                """,
                plot=self.kmer_signal_table(),
            )

            self.add_section(
                name="K-mer mode read parsing",
                anchor="qc3C-kmer-acceptance-plot",
                description="""
                    This figure displays a breakdown of proportion of parsed reads rejected due to various
                    criteria and the proportion that were accepted.
                """,
                plot=self.kmer_acceptance_plot(),
            )

            # TODO this figure seems of low value
            # self.add_section(
            #     name = 'K-mer mode junction proportion',
            #     anchor = 'qc3C-kmer-signal-plot',
            #     description = """
            #         This figure displays the proportion of inspected reads which contained a putative junction
            #         sequence, as compared against all reads inspected.
            #     """,
            #     plot = self.kmer_signal_plot()
            # )

            # TODO consider whether this figure should be displayed when there are non-trivial digests
            if self.do_digest_plot:
                self.add_section(
                    name="K-mer mode junction breakdown",
                    anchor="qc3C-kmer-junction-plot",
                    description="""
                        This figure displays the frequency at which a library's possible junction sequences
                        are actually observed in the reads. (_Trivial single-digests are ignored_)
                    """,
                    helptext="""
                        For trivial single-enzyme digests, there is only one possible junction sequence and so the
                        result for these experiments are not plotted. For dual-enzyme (such as Phase Genomics) there are
                        four potential junctions, while for dual-enzyme digests with one ambiguous site (such as Arima
                        Genomics) there are 16 possible junction sequences.

                        How efficiently the more complicated library protocols are at producing hybrid junctions is
                        possibly just a point of interest.

                        Junctions are named for which enzymes was responsible for creating the 5' and 3' ends.
                        E.g. `Sau3AI/MluCI` would involve two different enzymes, while `Sau3AI/Sau3AI` only one, as
                        would be the case in a single-enzyme digest. Proceeding the name is the actual junction
                        sequence.

                        The junctions are grouped by their 5' and then 3' enzyme, while the color spectrum used across
                        each bar aims to emphasise these enzymatic sources.

                    **Note:** in k-mer mode, the counts are not controlled for false positives.
                    """,
                    plot=self.kmer_junction_plot(),
                )

    @staticmethod
    def _drop_time(s):
        """remove time"""
        m = re.match(r"(.*) .+$", s)
        return s if m is None else m.group(1)

    @staticmethod
    def _drop_name(s):
        """drop leading program name from version string"""
        return s.split()[-1]

    def bam_runtime_table(self):
        config = {
            "id": "qc3C_bam_runtime_table",
            "namespace": "qc3C",
            "col1_header": "Sample",
            "scale": False,
            "title": "qc3C: BAM mode runtime parameters",
        }

        headers = {
            "b_run_timestamp": {
                "title": "Date",
                "description": "Analysis time stamp",
                "modify": MultiqcModule._drop_time,
                "hidden": True,
            },
            "b_mode": {"title": "Run Mode", "description": "Analysis mode used", "hidden": True},
            "b_min_mapq": {
                "title": "Min MapQ",
                "description": "Minimum accepted mapping quality",
                "min": 0,
                "format": "{:d}",
                "scale": False,
            },
            "b_enzymes": {"title": "Digest", "description": "Enzymes used in digest"},
            # 'b_seed': {
            #     'title': 'Seed',
            #     'description': 'Random seed',
            #     'format': '{:d}',
            #     'scale': False,
            #     'hidden': True,
            # },
            # 'b_max_obs': {
            #     'title': 'Max obs',
            #     'description': 'User specified maximum number of observations',
            #     'min': 0,
            #     'format': '{:,d}',
            #     'scale': 'OrRd',
            #     'modify': lambda x: 'n/a' if x == -1 else x
            # },
            "b_n_accepted_pairs": {
                "title": "Accepted pairs",
                "description": "Number of pairs accepted for analysis",
                "min": 0,
                "format": "{:,d}",
                "scale": "Oranges",
            },
            # 'b_sample_rate': {
            #     'title': 'Sample rate',
            #     'description': 'Sub-sampling probability',
            #     'min': 0,
            #     'max': 1,
            #     'format': '{:g}',
            #     'scale': 'Greys'
            # },
            # 'b_obs_insert_median': {
            #     'title': 'Insert median',
            #     'description': 'Estimated median insert size',
            #     'min': 0,
            #     'format': '{:,.0f}',
            #     'suffix': 'bp',
            #     'scale': 'Greens'
            # },
            "b_mean_readlen": {
                "title": "Read length",
                "description": "Average observed read length",
                "format": "{:,.0f}",
                "suffix": "bp",
                "scale": "Blues",
            },
            "b_obs_insert_mean": {
                "title": "Insert length",
                "description": "Inferred insert size",
                "min": 0,
                "format": "{:,.0f}",
                "suffix": "bp",
                "scale": "Purples",
            },
            "b_unobs_fraction": {
                "title": "Unobserved",
                "description": "Fraction of total insert extent that was unobservable",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Reds",
            },
            "b_p_read_thru": {
                "title": "Read-thru",
                "description": "Fraction of reads whose alignments end in a cutsite and whose sequence continues for the full junction",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
            },
        }

        return table.plot(self.qc3c_data["bam"], headers, config)

    def bam_longrange_plot(self):
        config = {
            "id": "qc3C_bam_longrange_plot",
            "title": "qc3C: BAM mode long-range pairs",
            "ylab": "Number of Reads",
            "cpswitch_counts_label": "Number of Reads",
        }

        categories = {
            "b_n_10k_": {"name": "10000 and further bp", "color": "#118AB2"},
            "b_n_5k_10k": {"name": "from 5000 to 10000 bp", "color": "#0CABAA"},
            "b_n_1k_5k": {"name": "from 1000 to 5000 bp", "color": "#07D1A1"},
            "b_n_short_inserts": {"name": "less than 1000 bp", "color": "#FFD166"},
            "b_n_trans_pairs": {"name": "trans pairs", "color": "#F78C6B"},
        }

        return bargraph.plot(self.qc3c_data["bam"], categories, config)

    def bam_acceptance_plot(self):
        config = {
            "id": "qc3C_bam_acceptance_plot",
            "title": "qc3C: BAM mode read parsing results",
            "ylab": "Number of Reads",
            "hide_empty": False,
            "cpswitch_counts_label": "Number of Reads",
        }

        categories = {
            "b_n_accepted_reads": {"name": "Accepted", "color": rev_8[-1]},
            "b_n_ref_term_reads": {"name": "Truncated", "color": rev_8[6]},
            "b_n_weak_mapping_reads": {"name": "Weak mapping", "color": rev_8[5]},
            "b_n_supplementary_reads": {"name": "Supplementary", "color": rev_8[4]},
            "b_n_secondary_reads": {"name": "Secondary", "color": rev_8[3]},
            "b_n_low_mapq_reads": {"name": "Low mapq", "color": rev_8[2]},
            "b_n_unmapped_reads": {"name": "Unmapped", "color": rev_8[1]},
            "b_n_skipped_reads": {"name": "Skipped", "color": rev_8[0]},
        }
        return bargraph.plot(self.qc3c_data["bam"], categories, config)

    def bam_signal_table(self):
        config = {"id": "qc3C_bam_signal_table", "namespace": "qc3C", "hide_empty": False, "col1_header": "Sample"}

        headers = {
            "b_unobs_fraction": {
                "title": "Unobservable extent",
                "description": "Estimated fraction of total fragment extent that was unobservable",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Reds",
            },
            "b_p_read_thru": {
                "title": "Read-thru",
                "description": "Fraction of reads whose alignments end in a cutsite and whose sequence continues for the full junction",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
            },
        }
        return table.plot(self.qc3c_data["bam"], headers, config)

    def bam_hicpro_table(self):
        config = {"id": "qc3C_bam_hicpro_table", "namespace": "qc3C", "hide_empty": False, "col1_header": "Sample"}

        headers = {
            "b_p_informative_fr": {"title": "Valid FR", "min": 0, "max": 100, "suffix": "%", "scale": "Greens"},
            "b_p_informative_rf": {"title": "Valid RF", "min": 0, "max": 100, "suffix": "%", "scale": "Greens"},
            "b_p_informative_ffrr": {
                "title": "Valid FF|RR",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
            },
            "b_p_uninformative_religation": {
                "title": "Religation",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Blues",
            },
            "b_p_uninformative_dangling_ends": {
                "title": "Dangling End",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Reds",
            },
            "b_p_uninformative_self_circle": {
                "title": "Self-circle",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Reds",
            },
            "b_p_uninformative_ffrr": {
                "title": "Invalid FF|RR",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Reds",
            },
        }
        return table.plot(self.qc3c_data["bam"], headers, config)

    def bam_valid_plot(self):
        config = {
            "id": "qc3C_bam_valid_plot",
            "title": "qc3C: BAM mode valid vs invalid HiC-Pro categories",
            "ylab": "Number of Reads",
            "hide_empty": False,
            "cpswitch_counts_label": "Number of Reads",
        }

        categories = {
            "b_n_informative_fr": {"name": "Valid FR", "color": "#41ab5d"},
            "b_n_informative_rf": {"name": "Valid RF", "color": "#74c476"},
            "b_n_informative_ffrr": {"name": "Valid FF|RR", "color": "#a1d99b"},
            "b_n_uninformative_religation": {"name": "Religation", "color": "#fcbba1"},
            "b_n_uninformative_dangling_ends": {"name": "Dangling End", "color": "#fc9272"},
            "b_n_uninformative_self_circle": {"name": "Self-circle", "color": "#fb6a4a"},
            "b_n_uninformative_ffrr": {"name": "Invalid FF|RR", "color": "#ef3b2c"},
        }
        return bargraph.plot(self.qc3c_data["bam"], categories, config)

    def bam_junction_plot(self):
        config = {
            "id": "qc3C_bam_junction_plot",
            "title": "qc3C: BAM mode read-thru ligation product frequency",
            "ylab": "Number of reads",
            "hide_empty": False,
            "use_legend": False,
        }

        categories = dict()
        for v in self.digest_junctions["bam"].values():
            for vi in v:
                categories[vi["name"]] = vi

        return bargraph.plot(self.qc3c_data["bam"], categories, config)

    def bam_fragment_histogram(self):
        median_lines = []
        for smpl in self.qc3c_data["bam"]:
            median_lines.append(
                {
                    "value": self.qc3c_data["bam"][smpl]["b_obs_insert_median"],
                    "color": "#D8E2DC",
                    "width": 2,
                    "dash": "dashdot",
                }
            )

        config = {
            "id": "qc3C_bam_fragment_histogram",
            "title": "qc3C: BAM mode distribution of pair separation",
            "cpswitch_counts_label": "Density",
            "logswitch": True,
            "logswitch_active": True,
            "logswitch_label": "Log10 [Density]",
            "xlog": True,
            "x_lines": median_lines,
            "xlab": "Log10 [Separation]",
            "ylab": "Density",
            "tt_label": "x:{point.x:.0f} bp, y:{point.y:.8f}",
        }

        data = {}
        for smpl in self.qc3c_data["bam"]:
            data[smpl] = self.qc3c_data["bam"][smpl]["frag_hist"]

        return linegraph.plot(data, config)

    def kmer_runtime_table(self):
        config = {
            "id": "qc3C_kmer_runtime_table",
            "namespace": "qc3C",
            "col1_header": "Sample",
            "scale": False,
            "title": "qc3C: K-mer mode runtime parameters",
        }

        headers = {
            "k_run_timestamp": {
                "title": "Date",
                "description": "Analysis time stamp",
                "modify": MultiqcModule._drop_time,
                "hidden": True,
            },
            "k_mode": {
                "title": "Run Mode",
                "description": "Analysis mode used",
                "hidden": True,
            },
            "k_kmer_size": {
                "title": "k",
                "description": "Library k-mer size",
                "min": 0,
                "format": "{:,d}",
                "scale": False,
                "hidden": True,
            },
            "k_enzymes": {"title": "Digest", "description": "Enzymes used in digest"},
            "k_n_accepted_reads": {
                "title": "Accepted reads",
                "description": "Number of reads accepted for analysis",
                "min": 0,
                "format": "{:,d}",
                "scale": "BuGn",
            },
            "k_mean_insert": {
                "title": "Insert length",
                "description": "User-specified insert size",
                "min": 0,
                "format": "{:,.0f}",
                "suffix": "bp",
                "scale": "Greens",
            },
            "k_mean_readlen": {
                "title": "Read length",
                "description": "Observed average read length",
                "format": "{:,.0f}",
                "suffix": "bp",
                "scale": "Blues",
            },
            "k_unobs_fraction": {
                "title": "Unobservable extent",
                "description": "Estimated mean of the unobservable portion of fragments",
                "shared_key": "unobs_mean",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Reds",
            },
        }
        return table.plot(self.qc3c_data["kmer"], headers, config)

    def kmer_signal_table(self):
        config = {
            "id": "qc3C_kmer_signal_table",
            "namespace": "qc3C",
            "col1_header": "Sample",
            "title": "qc3C: K-mer mode Hi-C fraction",
        }

        headers = {
            "k_raw_fraction": {
                "title": "Mean raw Hi-C fraction",
                "description": "Estimated mean of Hi-C fraction from only the observable extent",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Greens",
            },
            "k_adj_fraction": {
                "title": "Mean adjusted Hi-C fraction",
                "description": "Estimated mean of Hi-C fraction adjusted for unobservable extent",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "Blues",
            },
        }
        return table.plot(self.qc3c_data["kmer"], headers, config)

    def kmer_acceptance_plot(self):
        config = {
            "id": "qc3C_kmer_acceptance_plot",
            "title": "qc3C: K-mer mode read parsing results",
            "ylab": "Number of Reads",
            "hide_empty": False,
            "cpswitch_counts_label": "Number of Reads",
        }

        categories = {
            "k_n_accepted_reads": {"name": "Accepted", "color": rev_8[-1]},
            "k_n_low_cov": {"name": "Low cov", "color": rev_8[6]},
            "k_n_zero_cov": {"name": "Zero cov", "color": rev_8[5]},
            "k_n_high_cov": {"name": "High cov", "color": rev_8[4]},
            "k_n_ambiguous": {"name": "Ambiguous", "color": rev_8[3]},
            "k_n_no_flank": {"name": "No flank", "color": rev_8[2]},
            "k_n_too_short": {"name": "Too short", "color": rev_8[1]},
            "k_n_skipped": {"name": "Skipped", "color": rev_8[0]},
        }
        return bargraph.plot(self.qc3c_data["kmer"], categories, config)

    def kmer_signal_plot(self):
        config = {
            "id": "qc3C_kmer_signal_plot",
            "title": "qc3C: K-mer mode signal content",
            "ylab": "Number of Reads",
            "hide_empty": False,
            "stacking": None,
            "cpswitch": False,
            "cpswitch_c_active": False,
            "cpswitch_counts_label": "Number of Reads",
        }

        categories = {
            "k_raw_fraction": {"name": "Raw Hi-C fraction", "color": "#ef3b2c"},
            "k_adj_fraction": {"name": "Adjusted Hi-C fraction", "color": "#41ab5d"},
        }
        return bargraph.plot(self.qc3c_data["kmer"], categories, config)

    def kmer_junction_plot(self):
        config = {
            "id": "qc3C_kmer_frequency_plot",
            "title": "qc3C: K-mer mode putative ligation product frequency",
            "ylab": "Number of Reads",
            "hide_empty": False,
            "use_legend": False,
        }

        categories = dict()
        for _cat in self.digest_junctions["kmer"].values():
            for vi in _cat:
                categories[vi["name"]] = vi

        return bargraph.plot(self.qc3c_data["kmer"], categories, config)

    def parse_qc3c_log(self, f):
        def _none_to(x, y):
            return y if x is None else y

        parsed: Dict
        try:
            parsed_list = [json.loads(_l) for _l in f["f"]]
            if len(parsed_list) > 1:
                log.warning(
                    "Multiple records encountered in qc3C JSON file {}. Only the last report will be used".format(
                        os.path.join(f["root"], f["fn"])
                    )
                )
            parsed = parsed_list[-1]
        except json.JSONDecodeError:
            log.warning(f"Could not parse qc3C JSON: '{f['fn']}'")
            return

        s_name = self.clean_s_name(os.path.basename(os.path.abspath(f["root"])), f, root=os.path.dirname(f["root"]))
        if s_name in self.qc3c_data:
            log.debug(f"Duplicate sample name found! Overwriting: {f['s_name']}")

        # Add version info
        version = parsed["runtime_info"]["qc3C_version"].split()[-1]
        self.add_software_version(version, f["s_name"])

        try:
            analysis_mode = parsed["mode"]
            if "unobs_fraction" not in parsed and "unobs_frac" in parsed:
                parsed["unobs_fraction"] = parsed["unobs_frac"]

            if analysis_mode == "bam":
                if "vs_accepted" not in parsed["separation_bins"] and "zvs_accepted " in parsed["separation_bins"]:
                    parsed["separation_bins"]["vs_accepted"] = parsed["separation_bins"]["zvs_accepted"]

                # set some variables to shorten the lines below
                inf = parsed["classification"]["informative"]
                uninf = parsed["classification"]["uninformative"]
                n_cis_pairs = parsed["n_cis_pairs"]
                n_accepted_pairs = parsed["n_accepted_pairs"]
                n_paired_reads = n_accepted_pairs * 2

                self.qc3c_data["bam"][s_name] = {
                    "b_qc3C_version": parsed["runtime_info"]["qc3C_version"],
                    "b_run_timestamp": parsed["runtime_info"]["run_timestamp"],
                    "b_mode": parsed["mode"],
                    "b_enzymes": ", ".join(parsed["input_args"]["enzymes"]),
                    "b_seed": parsed["input_args"]["seed"],
                    "b_sample_rate": _none_to(parsed["input_args"]["sample_rate"], 1),
                    "b_max_obs": _none_to(parsed["input_args"]["max_obs"], -1),
                    "b_n_skipped_reads": parsed["n_skipped_reads"],
                    "b_n_unmapped_reads": parsed["n_unmapped"],
                    "b_n_analysed_reads": parsed["n_analysed_reads"],
                    "b_n_low_mapq_reads": parsed["n_low_mapq"],
                    "b_n_ref_len_reads": parsed["n_ref_len"],
                    "b_n_secondary_reads": parsed["n_secondary"],
                    "b_n_supplementary_reads": parsed["n_supplementary"],
                    "b_n_weak_mapping_reads": parsed["n_weak_mapping"],
                    "b_n_ref_term_reads": parsed["n_ref_term"],
                    "b_n_accepted_reads": parsed["n_accepted_reads"],
                    "b_obs_insert_mean": parsed["obs_insert_mean"],
                    "b_obs_insert_median": parsed["obs_insert_median"],
                    "b_mean_readlen": parsed["mean_readlen"],
                    "b_n_analysed_pairs": parsed["n_analysed_pairs"],
                    "b_n_accepted_pairs": parsed["n_accepted_pairs"],
                    "b_n_trans_pairs": parsed["n_trans_pairs"],
                    "b_p_trans_pairs": parsed["n_trans_pairs"] / n_accepted_pairs * 100,
                    "b_p_cis_pairs": n_cis_pairs / n_accepted_pairs * 100,
                    "b_unobs_fraction": parsed["unobs_fraction"] * 100,
                    "b_p_cs_start": parsed["digest_stats"]["cs_start"] / n_paired_reads * 100,
                    "b_p_cs_term": parsed["digest_stats"]["cs_term"] / n_paired_reads * 100,
                    "b_p_cs_full": parsed["digest_stats"]["cs_full"] / n_paired_reads * 100,
                    "b_p_read_thru": parsed["digest_stats"]["read_thru"] / n_paired_reads * 100,
                    "b_p_is_split": parsed["digest_stats"]["is_split"] / n_paired_reads * 100,
                    "b_adj_read_thru": parsed["digest_stats"]["read_thru"]
                    / n_paired_reads
                    * 100
                    * 1
                    / (1 - parsed["unobs_fraction"]),
                    "b_n_short_inserts": parsed["n_short_inserts"],
                    "b_p_short_inserts": parsed["n_short_inserts"] / n_accepted_pairs * 100,
                    "b_n_informative_fr": inf["fr"],
                    "b_n_informative_rf": inf["rf"],
                    "b_n_informative_ffrr": inf["ffrr"],
                    "b_n_uninformative_religation": uninf["religation"],
                    "b_n_uninformative_dangling_ends": uninf["dangling_ends"],
                    "b_n_uninformative_self_circle": uninf["self_circle"],
                    "b_n_uninformative_ffrr": uninf["ffrr"],
                    "b_p_informative_fr": inf["fr"] / n_cis_pairs * 100,
                    "b_p_informative_rf": inf["rf"] / n_cis_pairs * 100,
                    "b_p_informative_ffrr": inf["ffrr"] / n_cis_pairs * 100,
                    "b_p_uninformative_religation": uninf["religation"] / n_cis_pairs * 100,
                    "b_p_uninformative_dangling_ends": uninf["dangling_ends"] / n_cis_pairs * 100,
                    "b_p_uninformative_self_circle": uninf["self_circle"] / n_cis_pairs * 100,
                    "b_p_uninformative_ffrr": uninf["ffrr"] / n_cis_pairs * 100,
                    "b_n_1kb": parsed["separation_bins"]["counts"][0],
                    "b_n_5kb": parsed["separation_bins"]["counts"][1],
                    "b_n_10kb": parsed["separation_bins"]["counts"][2],
                    "b_p_1kb_vs_accepted": parsed["separation_bins"]["vs_accepted"][0],
                    "b_p_5kb_vs_accepted": parsed["separation_bins"]["vs_accepted"][1],
                    "b_p_10kb_vs_accepted": parsed["separation_bins"]["vs_accepted"][2],
                    "b_p_1kb_vs_cis": parsed["separation_bins"]["vs_all_cis"][0],
                    "b_p_5kb_vs_cis": parsed["separation_bins"]["vs_all_cis"][1],
                    "b_p_10kb_vs_cis": parsed["separation_bins"]["vs_all_cis"][2],
                    "b_n_1k_5k": parsed["separation_bins"]["counts"][0] - parsed["separation_bins"]["counts"][1],
                    "b_n_5k_10k": parsed["separation_bins"]["counts"][1] - parsed["separation_bins"]["counts"][2],
                    "b_n_10k_": parsed["separation_bins"]["counts"][2],
                }

                fhist = {}
                for x, y in itertools.zip_longest(
                    parsed["separation_histogram"]["mid_points"], parsed["separation_histogram"]["counts"]
                ):
                    fhist[float(x)] = float(y)
                self.qc3c_data["bam"][s_name]["frag_hist"] = fhist

            elif analysis_mode == "kmer":
                for k in "raw_fraction", "adj_fraction", "unobs_fraction":
                    if parsed[k] is None:
                        parsed[k] = "Error - adjusted value would exceed 100"
                    else:
                        parsed[k] = np.array(parsed[k]).mean() * 100

                self.qc3c_data["kmer"][s_name] = {
                    "k_qc3C_version": parsed["runtime_info"]["qc3C_version"],
                    "k_run_timestamp": parsed["runtime_info"]["run_timestamp"],
                    "k_mode": parsed["mode"],
                    "k_kmer_size": parsed["input_args"]["kmer_size"],
                    "k_enzymes": ", ".join(parsed["input_args"]["enzymes"]),
                    "k_seed": parsed["input_args"]["seed"],
                    "k_sample_rate": _none_to(parsed["input_args"]["sample_rate"], 1),
                    "k_max_freq": parsed["input_args"]["max_coverage"],
                    "k_mean_insert": parsed["input_args"]["mean_insert"],
                    "k_max_freq_quantile": parsed["input_args"]["max_freq_quantile"],
                    "k_max_obs": _none_to(parsed["input_args"]["max_obs"], -1),
                    "k_n_skipped": parsed["n_parsed_reads"] - parsed["n_analysed_reads"],
                    "k_n_analysed_reads": parsed["n_analysed_reads"],
                    "k_n_too_short": parsed["n_too_short"],
                    "k_n_no_flank": parsed["n_no_flank"],
                    "k_n_ambiguous": parsed["n_ambiguous"],
                    "k_n_high_cov": parsed["n_high_cov"],
                    "k_n_low_cov": parsed["n_low_cov"],
                    "k_n_zero_cov": parsed["n_zero_cov"],
                    "k_n_accepted_reads": parsed["n_accepted_reads"],
                    "k_n_without_junc": parsed["n_without_junc"],
                    "k_n_with_junc": parsed["n_with_junc"],
                    "k_mean_readlen": parsed["mean_readlen"],
                    "k_n_cs_start": parsed["cs_start"],
                    "k_raw_fraction": parsed["raw_fraction"],
                    "k_adj_fraction": parsed["adj_fraction"],
                    "k_unobs_fraction": parsed["unobs_fraction"],
                }

            # if any experiment contains a digest tht is non-trivial and produces
            # more than one possible junction sequence, prepare the supporting data
            # to render the junction frequency plot
            if len(parsed["junction_frequency"]) > 1:
                self.do_digest_plot = True
                log.debug(f"Enabled junction frequency plot for non-trivial digest: {f['root']}")
                # include the junction frequencies (1 or many depending on digest)
                self.qc3c_data[analysis_mode][s_name].update(parsed["junction_frequency"])

                # calculate the degeneracy of junction sequences per enzymatic combination (5p end =/= 3p end)
                # this can vary due to ambiguous bases in restriction site
                degen_count = {
                    f"{v['enz5p']}/{v['enz3p']}": 4 ** v["junction"].count("N")
                    for k, v in parsed["digestion"]["junctions"].items()
                }
                # get a palette for this series
                cols = color_picker(list(degen_count.values()))
                # sort these in accordance with that above
                juncs = np.sort(
                    np.array(
                        [(k.split(" ")[0], k) for k in parsed["junction_frequency"]],
                        dtype=np.dtype([("a", "S100"), ("b", "S100")]),
                    )
                )

                # keep a record of how these should be colored per sample
                self.digest_junctions[analysis_mode][s_name] = [
                    {"name": juncs[i][1].decode(), "color": cols[i]} for i in range(len(juncs))
                ]

            self.add_data_source(f, s_name, section=analysis_mode)

        except KeyError as ex:
            log.debug(
                "The entry {} was not found in the qc3C JSON file '{}', skipping sample {}".format(
                    str(ex), os.path.join(f["root"], f["fn"]), f["s_name"]
                )
            )
            return
