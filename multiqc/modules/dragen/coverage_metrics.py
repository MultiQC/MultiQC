# coding=utf-8


import itertools
import logging
import re
from collections import defaultdict
from typing import List

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import table
from multiqc.utils.util_functions import write_data_file

from .utils import Metric, make_headers

# Initialise the logger
log = logging.getLogger(__name__)

NAMESPACE = "DRAGEN coverage"

METRICS = [
    # id_in_data    title (display name)  gen_stats  cov_table  unit  description  precision
    # Read stats:
    Metric(
        "Average alignment coverage over {}",
        "Depth",
        "#",
        "#",
        "x",
        "Coverage depth over {}: number of uniquely mapped bases to {} divided by the number of sites in {}.",
    ),
    Metric(
        "Uniformity of coverage (PCT > 0.2*mean) over {}",
        "Uniformity (>0.2×mean)",
        "hid",
        "#",
        "%",
        "Percentage of sites with coverage greater than 20% of the mean coverage in region. Demonstrates the uniformity of coverage, the higher the better.",
    ),
    Metric(
        "Average chr X coverage over {}",
        "X cov",
        None,
        "hid",
        "x",
        "Average chromosome X coverage over {}. Calculated as the number of bases that aligned to the chromosome X (or to the intersection of chromosome X with the target region) divided by the total number of loci in the chromosome X (or the intersection with the target region). If there is no chromosome X in the reference genome or the region does not intersect chromosome X, this metric shows as NA.",
    ),
    Metric(
        "Average chr Y coverage over {}",
        "Y cov",
        None,
        "hid",
        "x",
        "Average chromosome Y coverage over {}. Calculated as the number of bases that aligned to the chromosome Y (or to the intersection of chromosome Y with the target region) divided by the total number of loci in the chromosome Y (or the intersection with the target region). If there is no chromosome Y in the reference genome or the region does not intersect chromosome Y, this metric shows as NA.",
    ),
    Metric(
        "Average mitochondrial coverage over {}",
        "MT cov",
        None,
        "hid",
        "x",
        "Average mitochondrial chromosome coverage over {}. Calculated as the number of bases that aligned to the mitochondrial chromosome (or to the intersection of mitochondrial chromosome with the target region) divided by the total number of loci in the mitochondrial chromosome (or the intersection with the target region). If there is no mitochondrial chromosome in the reference genome or the region does not intersect mitochondrial chromosome, this metric shows as NA.",
    ),
    Metric(
        "Average autosomal coverage over {}",
        "Mean aut cov",
        None,
        "hid",
        "x",
        "Average autosomal coverage over {}. Calculated as the number of bases that aligned to the autosomal loci in {} divided by the total number of loci in the autosomal loci in {}. If there is no autosomes in the reference genome, or the region does not intersect autosomes, this metric shows as NA.",
    ),
    Metric(
        "Median autosomal coverage over {}",
        "Med aut cov",
        "hid",
        "hid",
        "x",
        "Median alignment coverage over the autosomal loci in {}. If there is no autosome in the reference genome or the region does not intersect autosomes, this metric shows as NA.",
    ),
    Metric(
        "Mean/Median autosomal coverage ratio over {}",
        "Mean/med autosomal coverage",
        "hid",
        "#",
        "",
        "Mean autosomal coverage in {} divided by the median autosomal coverage in {}. If there is no autosome in the reference genome or the region does not intersect autosomes, this metric shows as NA.",
        precision=2,
    ),
    Metric(
        "PCT of {} with coverage [  1x: inf)",
        "⩾1x",
        "hid",
        "#",
        "%",
        "Percentage of sites in region with at least 1x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [  3x: inf)",
        "⩾3x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with at least 3x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 10x: inf)",
        "⩾10x",
        None,
        "#",
        "%",
        "Percentage of sites in region with at least 10x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 15x: inf)",
        "⩾15x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with at least 15x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 20x: inf)",
        "⩾20x",
        "#",
        "#",
        "%",
        "Percentage of sites in region with at least 20x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 50x: inf)",
        "⩾50x",
        "hid",
        "#",
        "%",
        "Percentage of sites in region with at least 50x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [100x: inf)",
        "⩾100x",
        "hid",
        "#",
        "%",
        "Percentage of sites in region with at least 100x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [  0x:  1x)",
        "0x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with no coverage.",
    ),
    Metric(
        "PCT of {} with coverage [  1x:  3x)",
        "1x..3x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with at least 1x but less than 3x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [  3x: 10x)",
        "3x..10x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with at least 3x but less than 10x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 10x: 15x)",
        "10x..15x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with at least 10x but less than 15x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 15x: 20x)",
        "15x..20x",
        None,
        "hid",
        "%",
        "Percentage of sites in region with at least 15x but less than 20x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 20x: 50x)",
        "20x..50x",
        "hid",
        "hid",
        "%",
        "Percentage of sites in region with at least 20x but less than 50x coverage.",
    ),
    Metric(
        "PCT of {} with coverage [ 50x:100x)",
        "50x..100x",
        "hid",
        "hid",
        "%",
        "Percentage of sites in region with at least 50x but less than 100x coverage.",
    ),
]

# For when target_bed or qc-coverage-region-x are available.
TARGET_COV_ONLY_METRICS = [
    Metric(
        "Aligned reads in {}",
        "Reads on target",
        "hid",
        "%",
        "reads",
        "Number of uniquely mapped reads to region relative to the number of uniquely mapped reads to the genome. When region is the target BED, this metric is equivalent to and replaces Capture Specificity based on target region.",
    ),
    Metric(
        "Aligned bases in {}",
        "Bases on target",
        "hid",
        "%",
        "bases",
        "Number of uniquely mapped bases to the region relative to the number of uniquely mapped bases to the genome.",
    ),
]

# Note only difference is 'Aligned reads in {}' in target cov but in wgs cov metrics its just 'Aligned reads'
WGS_COV_ONLY_METRICS = [
    Metric("Aligned reads", "Aln reads", "hid", "#", "reads", "Total number of aligned reads."),
    Metric("Aligned bases", "Aln bases", "hid", "#", "bases", "Total number of aligned bases."),
]

WGS_COV_METRICS = [
    Metric(
        m.id.replace("{}", "genome"),
        m.title.replace("{}", "genome"),
        in_genstats=m.in_genstats,
        in_own_tabl=m.in_own_tabl,
        unit=m.unit,
        descr=m.descr.replace("{}", "genome"),
        namespace=NAMESPACE,
        precision=m.precision,
    )
    for m in (WGS_COV_ONLY_METRICS + METRICS)
]

TARGET_COV_METRICS = [
    Metric(
        m.id.replace("{}", "target region"),
        m.title.replace("{}", "target region"),
        in_genstats=m.in_genstats,
        in_own_tabl=m.in_own_tabl,
        unit=m.unit,
        descr=m.descr.replace("{}", "target region"),
        namespace=NAMESPACE,
        precision=m.precision,
    )
    for m in (TARGET_COV_ONLY_METRICS + METRICS)
]


QC_REGION_COV_METRICS = [
    Metric(
        m.id.replace("{}", "QC coverage region"),
        m.title.replace("{}", "region"),
        in_genstats="hid",
        in_own_tabl=m.in_own_tabl,
        unit=m.unit,
        descr=m.descr.replace("{}", "region"),
        namespace=NAMESPACE,
        precision=m.precision,
    )
    for m in (TARGET_COV_ONLY_METRICS + METRICS)
]


class DragenCoverageMetrics(BaseMultiqcModule):
    def add_coverage_metrics_base_method(
        self,
        log_file_path: str,
        coverage_metrics_file_pattern: str,
        data_file_key: str,
        metrics: List,
        table_name: str,
        table_anchor: str,
    ):
        """
        :param log_file_path: One of
          'dragen/wgs_coverage_metrics',
          'dragen/qc_region_coverage_metrics'
          'dragen/target_bed_coverage_metrics'
        :param coverage_metrics_file_pattern: One of
          '(.*)\.wgs_coverage_metrics_?(tumor|normal)?.csv'
          '(.*)\.qc-coverage-region-([^_]+)_coverage_metrics_?(tumor|normal)?.csv'
          '(.*)\.target_bed_coverage_metrics_?(tumor|normal)?.csv'
        :param data_file_key: One of
          'dragen_wgs_cov_metrics'
          'dragen_qc_region_{region}_coverage_metrics'
          'dragen_target_bed_cov_metrics'
        :param metrics: One of
           WGS_COV_METRICS, QC_REGION_COV_METRICS, TARGET_COV_METRICS
        :param table_name: One of
          'WGS Coverage Metrics'
          'QC Region {region} Coverage Metrics'
          'Target Bed Coverage Metrics'
        :param table_anchor: One of
          'dragen-wgs-coverage-metrics'
          'dragen-qc-region-coverage-metrics-{region}'
          'dragen-target-bed-coverage-metrics'
        """
        data_by_phenotype_by_sample = defaultdict(dict)

        for f in self.find_log_files(log_file_path):
            sample_name, region_num, phenotype, data = parse_coverage_metrics(f, coverage_metrics_file_pattern)

            if sample_name is None:
                # Didn't match regex, continue
                continue

            # Region number is none if qc-region coverage metrics
            sample_name = self.clean_s_name(sample_name, f)
            if sample_name in data_by_phenotype_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {sample_name}")
            self.add_data_source(f, section="stats")
            data_by_phenotype_by_sample[sample_name].update({phenotype: data})

        # Filter to strip out ignored sample names:
        data_by_phenotype_by_sample = self.ignore_samples(data_by_phenotype_by_sample)

        # Merge tumor and normal data:
        data_by_sample = defaultdict(dict)
        for sn in data_by_phenotype_by_sample:
            for phenotype in data_by_phenotype_by_sample[sn]:
                new_sn = sn
                if phenotype == "normal":
                    new_sn = sn + "_normal"
                data_by_sample[new_sn] = data_by_phenotype_by_sample[sn][phenotype]

        if not data_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_sample, data_file_key)

        all_metric_names = set()
        for sn, sdata in data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)
        gen_stats_headers, own_tabl_headers = make_headers(all_metric_names, metrics)

        self.general_stats_addcols(data_by_sample, gen_stats_headers, namespace=NAMESPACE)

        return self._create_table(
            data_by_sample,
            own_tabl_headers,
            table_name,
            table_anchor,
        )

    def add_wgs_coverage_metrics(self):
        return self.add_coverage_metrics_base_method(
            log_file_path="dragen/wgs_coverage_metrics",
            coverage_metrics_file_pattern=r"(.*)\.wgs_coverage_metrics_?(tumor|normal)?.csv",
            data_file_key="dragen_wgs_cov_metrics",
            metrics=WGS_COV_METRICS,
            table_name="WGS Coverage Metrics",
            table_anchor="dragen-wgs-coverage-metrics",
        )

    def add_qc_region_coverage_metrics(self):
        data_by_sample_keys_set = set()

        for qc_coverage_num in range(1, 4):  # Only 3 qc coverage files available
            data_by_sample_keys_set |= self.add_coverage_metrics_base_method(
                log_file_path="dragen/qc_region_coverage_metrics",
                coverage_metrics_file_pattern=rf"(.*)\.qc-coverage-region-({qc_coverage_num})_coverage_metrics_?(tumor|normal)?.csv",
                data_file_key=f"dragen_qc_region_{qc_coverage_num}_coverage_metrics",
                metrics=QC_REGION_COV_METRICS,
                table_name=f"QC Region {qc_coverage_num} Coverage Metrics",
                table_anchor=f"dragen-qc-region-coverage-metrics-{qc_coverage_num}",
            )

        return data_by_sample_keys_set

    def add_target_bed_coverage_metrics(self):
        return self.add_coverage_metrics_base_method(
            log_file_path="dragen/target_bed_coverage_metrics",
            coverage_metrics_file_pattern=r"(.*)\.target_bed_coverage_metrics_?(tumor|normal)?.csv",
            data_file_key="dragen_target_bed_coverage_metrics",
            metrics=TARGET_COV_METRICS,
            table_name="Target Bed Coverage Metrics",
            table_anchor="dragen-target-bed-coverage-metrics",
        )

    def _create_table(self, data_by_sample, own_tabl_headers, table_name: str, table_anchor: str):
        self.add_section(
            name=table_name,
            anchor=table_anchor,
            description="""
            Coverage metrics over a region (where the region can be a target region,
            a QC coverage region, or the whole genome). Press the `Help` button for details.
            """,
            helptext="""
            The following criteria are used when calculating coverage:

            * Duplicate reads and clipped bases are ignored.
            * Only reads with `MAPQ` > `min MAPQ` and bases with `BQ` > `min BQ` are considered

            Considering only bases usable for variant calling, _i.e._ excluding:

            1. Clipped bases
            2. Bases in duplicate reads
            3. Reads with `MAPQ` < `min MAPQ` (default `20`)
            4. Bases with `BQ` < `min BQ` (default `10`)
            5. Reads with `MAPQ` = `0` (multimappers)
            6. Overlapping mates are double-counted
            """,
            plot=table.plot(data_by_sample, own_tabl_headers, pconfig={"namespace": NAMESPACE}),
        )
        return data_by_sample.keys()


def parse_coverage_metrics(f, file_regex):
    """
    T_SRR7890936_50pc.wgs_coverage_metrics_tumor.csv
    T_SRR7890936_50pc.wgs_coverage_metrics_normal.csv
    T_SRR7890936_50pc.wgs_coverage_metrics.csv
    T_SRR7890936_50pc.target_bed_coverage_metrics_tumor.csv
    T_SRR7890936_50pc.target_bed_coverage_metrics_normal.csv
    T_SRR7890936_50pc.target_bed_coverage_metrics.csv
    T_SRR7890936_50pc.qc-coverage-region-1_coverage_metrics_tumor.csv
    T_SRR7890936_50pc.qc-coverage-region-1_coverage_metrics_normal.csv
    T_SRR7890936_50pc.qc-coverage-region-1_coverage_metrics.csv


    The coverage metrics report outputs a _coverage_metrics.csv file, which provides metrics over a region,
    where the region can be the genome, a target region, or a QC coverage region. The first column of the output
    file contains the section name COVERAGE SUMMARY and the second column is empty for all metrics.

    The following criteria are used when calculating coverage:
    - Duplicate reads and clipped bases are ignored.
    - Only reads with MAPQ > min MAPQ and bases with BQ > min BQ are considered

    COVERAGE SUMMARY,,Aligned bases,250311219292
    COVERAGE SUMMARY,,Aligned bases in genome,250311219292,100.00
    COVERAGE SUMMARY,,Average alignment coverage over genome,82.09
    COVERAGE SUMMARY,,Uniformity of coverage (PCT > 0.2*mean) over genome,95.58
    COVERAGE SUMMARY,,PCT of genome with coverage [100x: inf),22.67
    COVERAGE SUMMARY,,PCT of genome with coverage [ 50x: inf),89.33
    COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: inf),95.44
    COVERAGE SUMMARY,,PCT of genome with coverage [ 15x: inf),95.68
    COVERAGE SUMMARY,,PCT of genome with coverage [ 10x: inf),95.92
    COVERAGE SUMMARY,,PCT of genome with coverage [  3x: inf),96.44
    COVERAGE SUMMARY,,PCT of genome with coverage [  1x: inf),96.91
    COVERAGE SUMMARY,,PCT of genome with coverage [  0x: inf),100.00
    COVERAGE SUMMARY,,PCT of genome with coverage [ 50x:100x),66.66
    COVERAGE SUMMARY,,PCT of genome with coverage [ 20x: 50x),6.10
    COVERAGE SUMMARY,,PCT of genome with coverage [ 15x: 20x),0.24
    COVERAGE SUMMARY,,PCT of genome with coverage [ 10x: 15x),0.24
    COVERAGE SUMMARY,,PCT of genome with coverage [  3x: 10x),0.52
    COVERAGE SUMMARY,,PCT of genome with coverage [  1x:  3x),0.47
    COVERAGE SUMMARY,,PCT of genome with coverage [  0x:  1x),3.09
    COVERAGE SUMMARY,,Average chr X coverage over genome,43.24
    COVERAGE SUMMARY,,Average chr Y coverage over genome,2.88
    COVERAGE SUMMARY,,Average mitochondrial coverage over genome,25675.54
    COVERAGE SUMMARY,,Average autosomal coverage over genome,85.00
    COVERAGE SUMMARY,,Median autosomal coverage over genome,83.79
    COVERAGE SUMMARY,,Mean/Median autosomal coverage ratio over genome,1.01
    COVERAGE SUMMARY,,Aligned reads,1689488809
    COVERAGE SUMMARY,,Aligned reads in genome,1689488809,100.00
    # where "genome" can be replaced with "target region" or "QC coverage region"

    Add metrics into gen stats table, plus add a dedicated table like VC
    """

    data = dict()

    for line in f["f"].splitlines():
        fields = line.split(",")
        metric = fields[2]
        value = fields[3]
        try:
            value = int(value)
        except ValueError:
            try:
                value = float(value)
            except ValueError:
                pass

        percentage = None
        if len(fields) > 4:  # percentage
            percentage = fields[4]
            try:
                percentage = float(percentage)
            except ValueError:
                pass

        data[metric] = value
        if percentage is not None:
            data[metric + " pct"] = percentage

    m = re.search(file_regex, f["fn"])

    if m is None:
        return None, None, None, None

    if len(m.groups()) == 3:
        sample, region_num, phenotype = m.group(1), m.group(2), m.group(3)
    elif len(m.groups()) == 2:
        sample, phenotype = m.group(1), m.group(2)
        region_num = None
    elif len(m.groups()) == 1:
        sample = m.group(1)
        region_num = None
        phenotype = None
    else:
        log.error(f"Could not get sample, phenotype from file name '{f['fn']}'")
        raise ValueError
    return sample, region_num, phenotype, data
