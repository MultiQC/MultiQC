import logging
from collections import defaultdict

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import table

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="PURPLE",
            anchor="purple",
            href="https://github.com/hartwigmedical/hmftools/",
            info="A purity, ploidy and copy number estimator for whole genome tumor data",
            extra="""
            PURPLE combines B-allele frequency (BAF), read depth ratios, somatic variants and
            structural variant breakpoints to estimate the purity and copy number profile
            of a tumor sample, and also predicts gender, the MSI status, tumor mutational
            load and burden, clonality and the whole genome duplication status.
            """,
            doi="10.1038/s41586-019-1689-y",
        )

        data_by_sample = defaultdict(dict)

        for f in self.find_log_files("purple/qc"):
            data = _parse_purple_qc(f)
            if data is not None:
                if f["s_name"] in data_by_sample:
                    log.debug(f"Duplicate PURPLE output prefix found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="stats")
                data_by_sample[f["s_name"]].update(data)

        for f in self.find_log_files("purple/purity"):
            data = _parse_purple_purity(f)
            if data is not None:
                if f["s_name"] in data_by_sample:
                    log.debug(f"Duplicate PURPLE output prefix found! Overwriting: {f['s_name']}")
                self.add_data_source(f, section="stats")
                data_by_sample[f["s_name"]].update(data)
                self.add_software_version(data.get("version"), f["s_name"])

        # Filter to strip out ignored sample names:
        data_by_sample = self.ignore_samples(data_by_sample)
        if not data_by_sample:
            raise ModuleNoSamplesFound
        log.info(f"Found {len(data_by_sample)} reports")

        # Write data to file
        self.write_data_file(data_by_sample, "purple")

        headers = _make_table_headers()

        self.add_section(
            name="PURPLE summary",
            anchor="purple-summary",
            description="""
            PURPLE summary. See details at the
            <a href=https://github.com/hartwigmedical/hmftools/tree/master/purity-ploidy-estimator#purity-file>
            documentation</a>.""",
            plot=table.plot(
                data_by_sample,
                headers,
                {
                    "id": "purple_summary",
                    "title": "PURPLE summary",
                },
            ),
            helptext="""
            1. **QC status**. Can fail for the following 3 reasons:
                * `FAIL_SEGMENT`: removed samples with more than 220 copy number segments unsupported
                    at either end by SV breakpoints. This step was added to remove samples with extreme GC bias,
                    with differences in depth of up to or in excess of 10x between high and low GC regions.
                    GC normalisation is unreliable when the corrections are so extreme so we filter.
                * `FAIL_GENDER`: if the AMBER and COBALT gender are inconsistent, we use the COBALT gender but
                    fail the sample.
                * `FAIL_DELETED_GENES`: we fail any sample with more than 280 deleted genes. This QC step was
                    added after observing that in a handful of samples with high MB scale positive GC bias we
                    sometimes systematically underestimate the copy number in high GC regions. This can lead us to
                    incorrectly infer homozygous loss of entire chromosomes, particularly on chromosomes 17 and 19.
            2. **Ploidy status**. Reflects how we have determined the purity of the sample:
                * `NORMAL`: could fix the purity using coverage and BAF alone
                * `HIGHLY_DIPLOID`: the fitted purity solution is highly diploid (> 95%) with a large
                    range of potential solutions, but somatic variants are unable to help either because they
                    were not supplied or because their implied purity was too low.
                * `SOMATIC`: somatic variants have improved the otherwise highly diploid solution
                * `NO_TUMOR`: PURPLE failed to find any aneuploidy and somatic variants were supplied but
                    there were fewer than 300 with observed VAF > 0.1.
            """,
        )

        gen_stat_cols = {k: v for k, v in headers.items() if k in ["ploidy", "purity"]}
        self.general_stats_addcols(data_by_sample, gen_stat_cols)


def _make_table_headers():
    headers = {
        "QCStatus": {
            "title": "QC Status",
            "description": "One of PASS, FAIL_SEGMENT, FAIL_GENDER, or FAIL_DELETED_GENES. "
            "For details, use the help button.",
            "scale": False,
        },
        "ploidy": {
            "title": "Ploidy",
            "description": "Average ploidy of the tumor sample after adjusting for purity",
            "scale": "RdYlGn",
            "min": 0,
        },
        "purity": {
            "title": "Purity",
            "description": "Purity of tumor in the sample",
            "scale": "RdYlGn",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "modify": lambda x: float(x) * 100.0,
        },
        "gender": {"title": "Gender", "description": "One of MALE, FEMALE or MALE_KLINEFELTER", "scale": False},
        "status": {
            "title": "Ploidy status",
            "description": "One of NORMAL, HIGHLY_DIPLOID, SOMATIC, or NO_TUMOR. For details, use the help button.",
            "scale": False,
        },
        "polyclonalProportion": {
            "title": "Polyclonal",
            "description": "Polyclonal proportion. Proportion of copy number regions that are more than 0.25 "
            "from a whole copy number",
            "scale": "RdYlGn",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "modify": lambda x: float(x) * 100.0,
        },
        "wholeGenomeDuplication": {
            "title": "WGD",
            "description": "Whole genome duplication. True if more than 10 autosomes have major allele ploidy > 1.5",
            "scale": False,
        },
        "msIndelsPerMb": {
            "title": "MS indel per Mb",
            "description": "Microsatellite indels per mega base",
            "scale": "RdYlGn",
            "hidden": True,
        },
        "msStatus": {
            "title": "MS status",
            "description": "Microsatellite status. One of MSI, MSS or UNKNOWN if somatic variants not supplied",
            "scale": False,
        },
        "tml": {
            "title": "TML",
            "description": "Tumor mutational load (# of missense variants in sample)",
            "scale": "RdYlGn",
            "hidden": True,
        },
        "tmlStatus": {
            "title": "TML status",
            "description": "Tumor mutational load status (# of missense variants in sample). One of HIGH, LOW or UNKNOWN "
            "if somatic variants not supplied",
            "scale": False,
        },
        "tmbPerMb": {
            "title": "TMB per Mb",
            "description": "Tumor mutational burden (# of passing variants) per mega base",
            "scale": "RdYlGn",
            "hidden": True,
        },
        "tmbStatus": {
            "title": "TMB status",
            "description": "Tumor mutational burden (# of passing variants per Mb) status. One of HIGH, LOW or UNKNOWN "
            "if somatic variants not supplied",
            "scale": False,
        },
    }
    return headers


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

    data = dict()
    for line in f["f"].splitlines():
        fields = line.strip().split("\t")
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

    header, values = [], []
    for i, line in enumerate(f["f"].splitlines()):
        fields = line.strip().split("\t")
        if i == 0:
            header = fields
        else:
            values = fields

    return dict(zip(header, values))
