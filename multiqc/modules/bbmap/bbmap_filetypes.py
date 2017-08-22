from collections import OrderedDict
from itertools import chain

from .plot_basic_hist import plot_basic_hist
from .plot_aqhist import plot_aqhist
from .plot_bhist import plot_bhist
#from .plot_bqhist import plot_bqhist
from .plot_idhist import plot_idhist
from .plot_mhist import plot_mhist
from .plot_qahist import plot_qahist
from .plot_qhist import plot_qhist

class slice2OrderedDict(object):
    def __getitem__(self, keys):
        return OrderedDict([(slice.start, slice.stop) for slice in keys])
odict = slice2OrderedDict()


file_types = {
    'aqhist': {
        'title': 'Average read quality',
        'descr': 'Histogram of average read quality.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'Quality':int,
            'count1':int,
            'fraction1':float,
            'count2':int,
            'fraction2':float
        ],
        'plot_func': plot_aqhist,
        'plot_params': {
            'yLog': True,
        }
    },
#basecov
#  - Coverage per base location.
#  - too big to interpret here
    'bhist' : {
        'title': 'Base composition',
        'descr': 'Base composition histogram by position.',
        'help_text': 'Relative composition',
        'cols': odict[
            'Pos':int, 'A':float, 'C':float, 'G':float, 'T':float, 'N':float
        ],
        'plot_func': plot_bhist,
        'plot_params': {}
    },
    'bincov': {
        'title': 'Binned coverage',
        'descr': 'Binned coverage per location (one line per X bases).',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean', 'STDev' ],
        'cols': odict[
            'RefName':str,
            'Cov':float,
            'Pos':int,
            'RunningPos':int
        ],
        'plot_func': plot_basic_hist,
        'plot_params': {},
        'not_implemented': ''
    },
    'bqhist': {
        'title': 'Base quality',
        'descr': 'Quality histogram designed for box plots.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'BaseNum':int,
            'count_1':int, 'min_1':int, 'max_1':int, 'mean_1':float,
            'Q1_1':int, 'med_1':int, 'Q3_1':int, 'LW_1':int, 'RW_1':int,
            'count_2':int, 'min_2':int, 'max_2':int, 'mean_2':float,
            'Q1_2':int, 'med_2':int, 'Q3_2':int, 'LW_2':int, 'RW_2':int
        ],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'covhist': {
        'title': 'Coverage histogram',
        'descr': 'Histogram of # occurrences of each depth level.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Coverage':int, 'numBases':int],
        'plot_func': plot_basic_hist,
        'plot_params': {
            'yLog': True,
        }
    },
    'covstats': {
        'title': 'Coverage stats',
        'descr': 'Per-scaffold coverage info.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'ID':str,
            'Avg_fold':float
        ],
        'extracols': odict[
            'Length':int,
            'Ref_GC':float,
            'Covered_percent':float,
            'Covered_bases':int,
            'Plus_reads':int,
            'Minus_reads':int,
            'Median_fold':int,
            'Read_GC':float,
            'Std_Dev':float
        ],
        'plot_func': plot_basic_hist,
        'plot_params': {},
        'not_implemented': ''
    },
    'ehist': {
        'title': 'Errors-per-read',
        'descr': 'Errors-per-read histogram.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Errors':int, 'Count':int ],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'gchist' : {
        'title': 'GC content',
        'descr': 'Read GC content histogram.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean', 'Median', 'Mode', 'STDev'],
        'cols': odict['GC':float, 'Count':int ],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'idhist': {
        'title': 'Identity histogram',
        'descr': 'Histogram of read count versus percent identity.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean_reads', 'Mean_bases',
                   'Median_reads', 'Median_bases',
                   'Mode_reads', 'Mode_bases',
                   'STDev_reads', 'STDev_bases' ],
        'cols': odict[
            'Identity':float, 'Reads':int, 'Bases':int
        ],
        'plot_func': plot_idhist,
        'plot_params': {}
    },
    'ihist': {
        'title': 'Insert sizes',
        'descr': 'Histogram of insert sizes (for paired reads).',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Mean', 'Median', 'STDev', 'PercentOfPairs'],
        'cols': odict['InsertSize':int, 'Count':int ],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'indelhist': {
        'title': 'Indel lengths',
        'descr': 'Indel length histogram.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Length':int, 'Deletions':int, 'Insertions':int],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'lhist' : {
        'title': 'Read lengths',
        'descr': 'Read length histogram.',
        'help_text': 'Placeholder help text.',
        'cols': odict['Length':int, 'Count':int ],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'mhist': {
        'title': 'Match, sub, del, and ins rates',
        'descr': 'Histogram of match, sub, del, and ins rates by read location.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'BaseNum':int,
            'Match1':float, 'Sub1':float, 'Del1':float, 'Ins1':float, 'N1':float, 'Other1':float,
            'Match2':float, 'Sub2':float, 'Del2':float, 'Ins2':float, 'N2':float, 'Other2':float
        ],
        'plot_func': plot_mhist,
        'plot_params': {}
    },
    'qahist': {
        'title': 'Quality accuracy',
        'descr': 'Quality accuracy histogram of error rates versus quality score.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['Deviation', 'DeviationSub' ],
        'cols': odict[
            'Quality':int, 'Match':int, 'Sub':int, 'Ins':int, 'Del':int,
            'TrueQuality':float, 'TrueQualitySub':float
        ],
        'plot_func': plot_qahist,
        'plot_params': {}
    },
    'qhist': {
        'title': 'Quality',
        'descr': 'Quality histogram by position.',
        'help_text': 'Placeholder help text.',
        'cols': odict[
            'BaseNum':int,
            'Read1_linear':float,
            'Read1_log':float,
            'Read1_measured':float,
            'Read2_linear':float,
            'Read2_log':float,
            'Read2_measured':float
        ],
        'plot_func': plot_qhist,
        'plot_params': {}
    },
    'rpkm': {
        'title': 'RPKM/FPKM',
        'descr': 'Per-scaffold RPKM/FPKM counts.',
        'help_text': 'Placeholder help text.',
        'kvrows': ['File', 'Reads', 'Mapped', 'RefSequences' ],
        'cols': odict[
            'Name':str,
            'Length':int, 'Bases':int, 'Coverage':float,
            'Reads':int, 'RPKM':float, 'Frags':int, 'FPKM':float
        ],
        'plot_func': plot_basic_hist,
        'plot_params': {},
        'not_implemented': '',
    },
    'statsfile_machine': {
        'title': 'General stats',
        'descr': 'General Stats',
        'help_text': 'Placeholder help text.',
        'cols': [],
        'plot_func': plot_basic_hist,
        'plot_params': {}
    },
    'statsfile': {
        'title': 'General stats',
        'descr': 'General Stats',
        'help_text': 'Placeholder help text.',
        'cols': [],
        'plot_params': {},
        'plot_func': plot_basic_hist,
        'not_implemented': ''
    }
}

statsfile_machine_keys = [
    'Reads_Used', 'Bases_Used', 'Reads/sec', 'kBases/sec',
    
    'R1_Mapped_Percent', 'R1_Unambiguous_Percent', 'R1_Mapped_Reads',
    'R1_Unambiguous_Reads', 'Mated_Pairs', 'Bad_Pairs', 'R1_Rescued',
    'Avg_Insert_Size', 'R1_Perfect_Best_Site', 'R1_Semiperfect_Site',
    'R1_Ambiguous_Mapping', 'R1_Low_Quality_Discards', 'R1_Match_Rate',
    'R1_Error_Rate', 'R1_Sub_Rate', 'R1_Del_Rate', 'R1_Ins_Rate',
    'R1_N_Rate', 'R1_Match_Count', 'R1_Error_Count', 'R1_Sub_Count',
    'R1_Del_Count', 'R1_Ins_Count', 'R1_N_Count',

    'R2_Mapped_Percent', 'R2_Unambiguous_Percent', 'R2_Mapped_Reads',
    'R2_Unambiguous_Reads', 'R2_Rescued', 'R2_Perfect_Best_Site',
    'R2_Semiperfect_Site', 'R2_Ambiguous_Mapping', 'R2_Low_Quality_Discards',
    'R2_Match_Rate', 'R2_Error_Rate', 'R2_Sub_Rate', 'R2_Del_Rate',
    'R2_Ins_Rate', 'R2_N_Rate', 'R2_Match_Count', 'R2_Error_Count',
    'R2_Sub_Count', 'R2_Del_Count', 'R2_Ins_Count', 'R2_N_Count',
    'R2_Mapped_Percent',
    
]

