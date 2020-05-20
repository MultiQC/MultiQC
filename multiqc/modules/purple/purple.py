#!/usr/bin/env python

""" MultiQC module to parse QC output from PURPLE """
import re
from collections import OrderedDict, defaultdict
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table

# Initialise the logger
import logging
log = logging.getLogger(__name__)

class MultiqcModule(BaseMultiqcModule):
    """
    PURPLE is a purity ploidy estimator. It combines B-allele frequency (BAF) from AMBER,
    read depth ratios from COBALT, somatic variants and structural variants to estimate the
    purity and copy number profile of a tumor sample.
    """

    def __init__(self):

        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name='PURPLE',
            anchor='purple',
            href="https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator",
            info="""is a purity ploidy estimator. It combines B-allele frequency (BAF) from AMBER,
                    read depth ratios from COBALT, somatic variants and structural variants to estimate the
                    purity and copy number profile of a tumor sample, as well as the MSI and the TMB status."""
        )

        data_by_sample = defaultdict(dict)

        for f in self.find_log_files('purple/qc'):
            data = _parse_purple_qc(f)
            if data is not None:
                if f['s_name'] in data_by_sample:
                    log.debug('Duplicate PURPLE output prefix found! Overwriting: {}'.format(f['s_name']))
                self.add_data_source(f, section='stats')
                data_by_sample[f['s_name']].update(data)

        for f in self.find_log_files('purple/purity'):
            data = _parse_purple_purity(f)
            if data is not None:
                if f['s_name'] in data_by_sample:
                    log.debug('Duplicate PURPLE output prefix found! Overwriting: {}'.format(f['s_name']))
                self.add_data_source(f, section='stats')
                data_by_sample[f['s_name']].update(data)

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            raise UserWarning
        log.info("Found {} reports".format(len(data_by_sample)))

        self._general_stats_table(data_by_sample)
        self._purple_stats_table(data_by_sample)

    def _general_stats_table(self, data_by_sample):
        headers = OrderedDict()
        headers['QCStatus'] = {
            'title': 'PURPLE QC',
            'description': 'PURPLE QC status (PASS, FAIL_SEGMENT, FAIL_GENDER, FAIL_DELETED_GENES).',
            'scale': False
        }
        self.general_stats_addcols(data_by_sample, headers)

    def _purple_stats_table(self, data_by_sample):
        headers = OrderedDict()
        headers['QCStatus'] = {
            'title': 'PURPLE QC',
            'description': 'PURPLE QC status (PASS, FAIL_SEGMENT, FAIL_GENDER, FAIL_DELETED_GENES).',
            'scale': False
        }
        headers['ploidy'] = {
            'title': 'Ploidy',
            'description': 'Average ploidy of the tumor sample after adjusting for purity',
            'scale': 'RdYlGn',
            'min': 0,
        }
        headers['purity'] = {
            'title': 'Purity',
            'description': 'Purity of tumor in the sample',
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'modify': lambda x: float(x) * 100.0,
        }
        headers['gender'] = {
            'title': 'Gender',
            'description': 'One of MALE, FEMALE or MALE_KLINEFELTER',
            'scale': False
        }
        headers['status'] = {
            'title': 'Ploidy status',
            'description': 'One of NORMAL, HIGHLY_DIPLOID, SOMATIC or NO_TUMOR',
            'scale': False
        }
        headers['polyclonalProportion'] = {
            'title': 'Polyclonal',
            'description': 'Polyclonal proportion. Proportion of copy number regions that are more than 0.25 from a whole copy number',
            'scale': 'RdYlGn',
            'min': 0,
            'max': 100,
            'suffix': '%',
            'modify': lambda x: float(x) * 100.0,
        }
        headers['wholeGenomeDuplication'] = {
            'title': 'WGD',
            'description': 'Whole genome duplication. True if more than 10 autosomes have major allele ploidy > 1.5',
            'scale': False
        }
        headers['msIndelsPerMb'] = {
            'title': 'MS indel per Mb',
            'description': 'Microsatellite indels per mega base',
            'scale': 'RdYlGn',
            'hidden': True,
        }
        headers['msStatus'] = {
            'title': 'MS status',
            'description': 'Microsatellite status. One of MSI, MSS or UNKNOWN if somatic variants not supplied',
            'scale': False
        }
        headers['tml'] = {
            'title': 'TML',
            'description': 'Tumor mutational load',
            'scale': 'RdYlGn',
            'hidden': True,
        }
        headers['tmlStatus'] = {
            'title': 'TML status',
            'description': 'Tumor mutational load status. One of HIGH, LOW or UNKNOWN if somatic variants not supplied',
            'scale': False
        }
        headers['tmbPerMb'] = {
            'title': 'TMB per Mb',
            'description': 'Tumor mutational burden per mega base',
            'scale': 'RdYlGn',
            'hidden': True,
        }
        headers['tmbStatus'] = {
            'title': 'TMB status',
            'description': 'Tumor mutational burden status. One of HIGH, LOW or UNKNOWN if somatic variants not supplied',
            'scale': False
        }
        self.add_section (
            name='PURPLE summary',
            anchor='purple-summary',
            description="PURPLE summary. See details at the "
                        "<a href=https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#purity-file>"
                        "documentation</a>.",
            plot=table.plot(data_by_sample, headers, {
                'id': 'purple_summary',
                'namespace': 'PURPLE',
                'title': 'PURPLE summary',
            })
        )


def _parse_purple_qc(f):
    """
    $ cat <sample>.purple.qc
    QCStatus        FAIL_DELETED_GENES
    SegmentPass     true
    GenderPass      true
    DeletedGenesPass        false
    SegmentScore    0
    UnsupportedSegments     0
    Ploidy  2.0036
    AmberGender     MALE
    CobaltGender    MALE
    DeletedGenes    5529
    """

    f['s_name'] = re.search(r'(.*).purple.qc', f['fn']).group(1)

    data = dict()
    for line in f['f'].splitlines():
        fields = line.strip().split('\t')
        if len(fields) == 2:
            data[fields[0]] = fields[1]
    return data


def _parse_purple_purity(f):
    """
    $ cat <sample>.purple.purity.tsv
    purity  normFactor  score   diploidProportion  ploidy  gender  status  polyclonalProportion  minPurity  maxPurity \
    minPloidy  maxPloidy  minDiploidProportion  maxDiploidProportion  version  somaticPenalty  wholeGenomeDuplication \
    msIndelsPerMb         msStatus  tml  tmlStatus  tmbPerMb            tmbStatus
    0.6300  1.1600      0.5126  0.0000             2.0036  MALE    NORMAL  0.0000                0.6000     0.6400 \
    1.9480     2.0037     0.0000                0.0000                2.40     0.0000          false \
    0.012941587967820916  MSS       0.0  LOW        1.0514165792235046  LOW
    """

    f['s_name'] = re.search(r'(.*).purple.purity.tsv', f['fn']).group(1)

    header, values = [], []
    for i, line in enumerate(f['f'].splitlines()):
        fields = line.strip().split('\t')
        if i == 0:
            header = fields
        else:
            values = fields

    return dict(zip(header, values))






