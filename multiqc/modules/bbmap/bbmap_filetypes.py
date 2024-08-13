from typing import Dict

from multiqc import config

from .plot_aqhist import plot_aqhist
from .plot_basic_hist import plot_basic_hist
from .plot_bhist import plot_bhist
from .plot_bqhist import plot_bqhist
from .plot_covhist import plot_covhist
from .plot_idhist import plot_idhist
from .plot_ihist import plot_ihist
from .plot_indelhist import plot_indelhist
from .plot_mhist import plot_mhist
from .plot_qahist import plot_qahist
from .plot_qchist import plot_qchist
from .plot_qhist import plot_qhist


section_order = [
    "stats",
    "covhist",
    "covstats",
    "bincov",
    "bqhist",
    "idhist",
    "indelhist",
    "mhist",
    "qahist",
    "qchist",
    "qhist",
    "aqhist",
    "ehist",
    "lhist",
    "ihist",
    "gchist",
    "bhist",
    "rpkm",
    "statsfile_machine",
    "statsfile",
]
file_types: Dict = {
    "stats": {
        "title": "BBDuk filtering statistics",
        "descr": "Proportion of reads that matched adapters/contaminants.",
        "help_text": "",
        "kvrows": ["Total", "Matched"],
        "kv_descriptions": {
            "Total": (
                "Total number of reads processed",
                {
                    "description": f"Aligned Reads ({config.read_count_desc})",
                    "shared_key": "read_count",
                    "modify": lambda x: x * config.read_count_multiplier,
                    "scale": "PuBu",
                    "hidden": True,
                },
            ),
            "Matched": (
                "Total number of reads matching adapters/contaminants",
                {
                    "description": f"Aligned Reads ({config.read_count_desc})",
                    "shared_key": "read_count",
                    "modify": lambda x: x * config.read_count_multiplier,
                    "scale": "Reds",
                    "hidden": True,
                },
            ),
            "Percent filtered": (
                "Proportion of reads filtered, matching adapters/contaminants",
                {"max": 100, "min": 0, "scale": "OrRd", "suffix": "%", "hidden": False},
            ),
        },
        "cols": {
            "Name": str,
            "Reads": int,
            "ReadsPct": lambda v: float(v.strip("%")),
        },
        "extracols": {
            "Bases": int,
            "BasesPct": float,
        },
        "plot_func": None,  # Plotting for 'stats' not implemented
        "plot_params": {},
    },
    "aqhist": {
        "title": "Read quality",
        "descr": "Histogram of average read qualities (`aqhist`). "
        "Plot shows the number of reads at each quality score.",
        "help_text": "",
        "cols": {"Quality": int, "count1": int, "fraction1": float, "count2": int, "fraction2": float},
        "plot_func": plot_aqhist,
        "plot_params": {
            "x_bands": [
                {"from": 28, "to": 100, "color": "#c3e6c3"},
                {"from": 20, "to": 28, "color": "#e6dcc3"},
                {"from": 0, "to": 20, "color": "#e6c3c3"},
            ],
            "ylog": True,
        },
    },
    # basecov
    #  - Coverage per base location.
    #  - too big to interpret here
    "bhist": {
        "title": "Base composition",
        "descr": "Base composition histogram by position (`bhist`). "
        "The plot shows the percentage of `G+C`, `A+T`, and `N` bases "
        "for each position in the reads.",
        "help_text": "Relative composition",
        "cols": {"Pos": int, "A": float, "C": float, "G": float, "T": float, "N": float},
        "plot_func": plot_bhist,
        "plot_params": {},
    },
    "bincov": {
        "title": "Binned coverage",
        "descr": "Binned coverage per location, one line per X bases (`bincov`).",
        "help_text": "",
        "kvrows": ["Mean", "STDev"],
        "cols": {"RefName": str, "Cov": float, "Pos": int, "RunningPos": int},
        "plot_func": plot_basic_hist,
        "plot_params": {},
        "not_implemented": "",
    },
    "bqhist": {
        "title": "Base quality",
        "descr": "Quality histogram designed for box plots (`bqhist`). "
        "Refer to original source files for complete boxplot data. "
        "Plot shows mean base quality for each read position. ",
        "help_text": "",
        "cols": {
            "BaseNum": int,
            "count_1": int,
            "min_1": int,
            "max_1": int,
            "mean_1": float,
            "Q1_1": int,
            "med_1": int,
            "Q3_1": int,
            "LW_1": int,
            "RW_1": int,
            "count_2": int,
            "min_2": int,
            "max_2": int,
            "mean_2": float,
            "Q1_2": int,
            "med_2": int,
            "Q3_2": int,
            "LW_2": int,
            "RW_2": int,
        },
        "plot_func": plot_bqhist,
        "plot_params": {},
    },
    "covhist": {
        "title": "Coverage histogram",
        "descr": "Histogram of number of occurrences of each coverage depth level (`covhist`). "
        "Note that lines have been smoothed to 400 points; "
        "higher resolution data might be available in the original data source. ",
        "help_text": "",
        "cols": {"Coverage": int, "numBases": int},
        "plot_func": plot_covhist,
        "plot_params": {
            "ylog": True,
        },
    },
    "covstats": {
        "title": "Coverage stats",
        "descr": "Per-scaffold coverage info (`covstats`).",
        "help_text": "",
        "cols": {"ID": str, "Avg_fold": float},
        "extracols": {
            "Length": int,
            "Ref_GC": float,
            "Covered_percent": float,
            "Covered_bases": int,
            "Plus_reads": int,
            "Minus_reads": int,
            "Median_fold": int,
            "Read_GC": float,
            "Std_Dev": float,
        },
        "plot_func": plot_basic_hist,
        "plot_params": {},
        "not_implemented": "",
    },
    "ehist": {
        "title": "Errors-per-read",
        "descr": "Errors-per-read histogram (`ehist`). ",
        "help_text": "",
        "cols": {"Errors": int, "Count": int},
        "plot_func": plot_basic_hist,
        "plot_params": {"xlab": "Errors", "ylab": "# Reads"},
    },
    "gchist": {
        "title": "GC content",
        "descr": "Read GC content histogram (`gchist`).",
        "help_text": "",
        "kvrows": ["Mean", "Median", "Mode", "STDev"],
        "kv_descriptions": {
            "Mean": ("Average GC content", {}),
            "Median": ("Median GC content", {}),
            "Mode": ("The most commonly occuring value of the GC content distribution", {}),
            "STDev": ("Standard deviation of average GC content", {}),
        },
        "cols": {"GC": float, "Count": int},
        "plot_func": plot_basic_hist,
        "plot_params": {
            "xlab": "Proportion GC",
            "ylab": "# Reads",
            "xsuffix": "%",
            "xmin": 0,
            "xmax": 100,
        },
    },
    "idhist": {
        "title": "Identity histogram",
        "descr": "Histogram of read count versus percent base pair identity " "of aligned reads (`idhist`).",
        "help_text": "",
        "kvrows": [
            "Mean_reads",
            "Mean_bases",
            "Median_reads",
            "Median_bases",
            "Mode_reads",
            "Mode_bases",
            "STDev_reads",
            "STDev_bases",
        ],
        "kv_descriptions": {
            "Mean_reads": ("Average percent identity of aligned reads", {}),
            "Mean_bases": ("Average percent identity of aligned bases", {}),
            "Median_reads": ("Median percent identity of aligned reads", {}),
            "Median_bases": ("Median percent identity of aligned bases", {}),
            "Mode_reads": (
                "The most commonly occuring average value amongst aligned reads (i.e. the mode of the distribution)",
                {},
            ),
            "Mode_bases": (
                "The most commonly occuring average value amongst aligned bases (i.e. the mode of the distribution)",
                {},
            ),
            "STDev_reads": (
                "The standard deviation of the average percent identity distribution for aligned reads",
                {},
            ),
            "STDev_bases": (
                "The standard deviation of the average percent identity distribution for aligned bases",
                {},
            ),
        },
        "cols": {"Identity": float, "Reads": int, "Bases": int},
        "plot_func": plot_idhist,
        "plot_params": {},
    },
    "ihist": {
        "title": "Insert sizes",
        "descr": "Histogram of computed insert sizes, for paired reads (`ihist`). "
        "Plotted data has been cut off at 99% to prevent long tails; "
        "Complete data available in original source files.",
        "help_text": "The insert size is the length of the sequence between the "
        "sequencing adapters, which for most common insert sizes is "
        "longer than the sum of both read pairs. "
        "In some cases, the insert size is shorter than the length "
        "of the two read pairs combined, resulting in an insert size "
        "shorter than the sum of the length of the reads pairs.",
        "kvrows": ["Mean", "Median", "STDev", "PercentOfPairs"],
        "kv_descriptions": {
            "Mean": ("Average insert length", {}),
            "Median": ("Median insert length", {}),
            "STDev": ("Standard deviation of insert size length distribution", {}),
            "PercentOfPairs": ("", {}),
        },
        "cols": {"InsertSize": int, "Count": int},
        "plot_func": plot_ihist,
        "plot_params": {},
    },
    "indelhist": {
        "title": "Indel lengths",
        "descr": "Indel length histogram (`indelhist`). "
        "The plots show the number of observed insertions and deletions, "
        "for each insertion and deletion length.",
        "help_text": "",
        "cols": {"Length": int, "Deletions": int, "Insertions": int},
        "plot_func": plot_indelhist,
        "plot_params": {},
    },
    "lhist": {
        "title": "Read lengths",
        "descr": "Read length histogram (`lhist`).",
        "help_text": "",
        "cols": {"Length": int, "Count": int},
        "plot_func": plot_basic_hist,
        "plot_params": {"xlab": "Read length (base pairs)", "xsuffix": "bp", "ylab": "# Reads"},
    },
    "mhist": {
        "title": "Match, substitution, deletion, and insertion rates",
        "descr": "Histogram of match, substitution, deleletion, " "and insertion rates by read location (`mhist`).",
        "help_text": "",
        "cols": {
            "BaseNum": int,
            "Match1": float,
            "Sub1": float,
            "Del1": float,
            "Ins1": float,
            "N1": float,
            "Other1": float,
            "Match2": float,
            "Sub2": float,
            "Del2": float,
            "Ins2": float,
            "N2": float,
            "Other2": float,
        },
        "plot_func": plot_mhist,
        "plot_params": {},
    },
    "qahist": {
        "title": "Quality accuracy",
        "descr": "Base quality accuracy histogram of error rates versus quality score (`qahist`). "
        "The plots show the observed count of each type of alignment by base quality score.",
        "help_text": "",
        "kvrows": ["Deviation", "DeviationSub"],
        "kv_descriptions": {
            "Deviation": ("", {}),
            "DeviationSub": ("", {}),
        },
        "cols": {
            "Quality": int,
            "Match": int,
            "Sub": int,
            "Ins": int,
            "Del": int,
            "TrueQuality": float,
            "TrueQualitySub": float,
        },
        "plot_func": plot_qahist,
        "plot_params": {},
    },
    "qchist": {
        "title": "Count of bases with each quality value",
        "descr": "Histogram of base qualities (`qchist`). "
        "Plot shows the number of bases at each quality score. Zero counts are shown as `0.1` due to log axis.",
        "help_text": "",
        "cols": {"Quality": int, "count1": int, "fraction1": float},
        "plot_func": plot_qchist,
        "plot_params": {
            "x_bands": [
                {"from": 30, "to": 100, "color": "#c3e6c3"},
                {"from": 20, "to": 30, "color": "#e6dcc3"},
                {"from": 0, "to": 20, "color": "#e6c3c3"},
            ],
            "ylog": True,
            "xlab": "Phred Score",
            "ylab": "Counts",
        },
    },
    "qhist": {
        "title": "Sequence Quality Histograms",
        "descr": "Quality histogram by position (`qhist`). "
        "The plots show the average quality for each position in the reads, "
        "using the linear values, logarithmically scaled values, and the "
        "actual measured qualities based on the alignments.",
        "help_text": "",
        "cols": {
            "BaseNum": int,
            "Read1_linear": float,
            "Read1_log": float,
            "Read1_measured": float,
            "Read2_linear": float,
            "Read2_log": float,
            "Read2_measured": float,
        },
        "plot_func": plot_qhist,
        "plot_params": {},
    },
    "rpkm": {
        "title": "RPKM/FPKM",
        "descr": "Per-scaffold RPKM/FPKM counts (`rpkm`).",
        "help_text": "",
        "kvrows": ["File", "Reads", "Mapped", "RefSequences"],
        "kv_descriptions": {
            "File": ("", {}),
            "Reads": ("", {}),
            "Mapped": ("", {}),
            "RefSequences": ("", {}),
        },
        "cols": {
            "Name": str,
            "Length": int,
            "Bases": int,
            "Coverage": float,
            "Reads": int,
            "RPKM": float,
            "Frags": int,
            "FPKM": float,
        },
        "plot_func": plot_basic_hist,
        "plot_params": {},
        "not_implemented": "",
    },
    "statsfile_machine": {
        "title": "General stats",
        "descr": "General Stats (`statsfile_machine`).",
        "help_text": "",
        "cols": [],
        "plot_func": plot_basic_hist,
        "plot_params": {},
    },
    "statsfile": {
        "title": "General stats",
        "descr": "General Stats (`statsfile`).",
        "help_text": "",
        "cols": [],
        "plot_params": {},
        "plot_func": plot_basic_hist,
        "not_implemented": "",
    },
}

statsfile_machine_keys = [
    "Reads_Used",
    "Bases_Used",
    "Reads/sec",
    "kBases/sec",
    "R1_Mapped_Percent",
    "R1_Unambiguous_Percent",
    "R1_Mapped_Reads",
    "R1_Unambiguous_Reads",
    "Mated_Pairs",
    "Bad_Pairs",
    "R1_Rescued",
    "Avg_Insert_Size",
    "R1_Perfect_Best_Site",
    "R1_Semiperfect_Site",
    "R1_Ambiguous_Mapping",
    "R1_Low_Quality_Discards",
    "R1_Match_Rate",
    "R1_Error_Rate",
    "R1_Sub_Rate",
    "R1_Del_Rate",
    "R1_Ins_Rate",
    "R1_N_Rate",
    "R1_Match_Count",
    "R1_Error_Count",
    "R1_Sub_Count",
    "R1_Del_Count",
    "R1_Ins_Count",
    "R1_N_Count",
    "R2_Mapped_Percent",
    "R2_Unambiguous_Percent",
    "R2_Mapped_Reads",
    "R2_Unambiguous_Reads",
    "R2_Rescued",
    "R2_Perfect_Best_Site",
    "R2_Semiperfect_Site",
    "R2_Ambiguous_Mapping",
    "R2_Low_Quality_Discards",
    "R2_Match_Rate",
    "R2_Error_Rate",
    "R2_Sub_Rate",
    "R2_Del_Rate",
    "R2_Ins_Rate",
    "R2_N_Rate",
    "R2_Match_Count",
    "R2_Error_Count",
    "R2_Sub_Count",
    "R2_Del_Count",
    "R2_Ins_Count",
    "R2_N_Count",
    "R2_Mapped_Percent",
]
