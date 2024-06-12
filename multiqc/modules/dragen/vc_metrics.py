import logging

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

from .utils import Metric, exist_and_number, make_headers

# Initialise the logger
log = logging.getLogger(__name__)


NAMESPACE = "Variant calling"


class DragenVCMetrics(BaseMultiqcModule):
    def add_vc_metrics(self):
        vc_data_by_sample = dict()
        for f in self.find_log_files("dragen/vc_metrics"):
            data = parse_vc_metrics_file(f)
            s_name = f["s_name"]
            if s_name in vc_data_by_sample:
                log.debug(f"Duplicate DRAGEN vc_metrics file was found! Overwriting sample {s_name}")
            self.add_data_source(f, s_name=s_name, section="vc_metrics")
            vc_data_by_sample[s_name] = data

        gvcf_data_by_sample = dict()
        for f in self.find_log_files("dragen/gvcf_metrics"):
            data = parse_vc_metrics_file(f)
            s_name = f["s_name"]
            if s_name in gvcf_data_by_sample:
                log.debug(f"Duplicate DRAGEN gvcf_metrics file was found! Overwriting sample {s_name}")
            self.add_data_source(f, s_name=s_name, section="gvcf_metrics")
            gvcf_data_by_sample[s_name] = data

        # Filter to strip out ignored sample names:
        vc_data_by_sample = self.ignore_samples(vc_data_by_sample)
        gvcf_data_by_sample = self.ignore_samples(gvcf_data_by_sample)
        if not vc_data_by_sample:
            return set()

        # Write data to file
        self.write_data_file(vc_data_by_sample, "dragen_vc_metrics")
        if gvcf_data_by_sample:
            self.write_data_file(gvcf_data_by_sample, "dragen_gvcf_metrics")

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        all_metric_names = set()
        for sn, sdata in vc_data_by_sample.items():
            for m in sdata.keys():
                all_metric_names.add(m)

        gen_stats_headers, table_headers = make_headers(all_metric_names, VC_METRICS)

        self.general_stats_addcols(vc_data_by_sample, gen_stats_headers, namespace=NAMESPACE)

        self.add_section(
            name="Variant calling",
            anchor="dragen-vc-metrics",
            description="""
            Metrics are reported for each sample in multi-sample VCF
            files. Based on the run case, metrics are reported either as standard
            VARIANT CALLER or JOINT CALLER. All metrics are reported for post-filter VCFs,
            except for the "Filtered" metrics which represent how many variants were filtered out
            from pre-filter VCF to generate the post-filter VCF.
            """,
            plot=table.plot(
                vc_data_by_sample,
                table_headers,
                pconfig={
                    "id": "dragen-vc-metrics-table",
                    "namespace": NAMESPACE,
                    "title": "DRAGEN: Variant calling metrics",
                },
            ),
        )

        if gvcf_data_by_sample:
            self.add_section(
                name="GVCF metrics",
                anchor="dragen-gvcf-metrics",
                description="""
                Metrics are calculated for each sample in a corresponding GVCF file.
                """,
                plot=table.plot(
                    gvcf_data_by_sample,
                    table_headers,
                    pconfig={
                        "id": "dragen-gvcf-metrics-table",
                        "namespace": NAMESPACE,
                        "title": "DRAGEN: GVCF metrics",
                    },
                ),
            )

        return vc_data_by_sample.keys()


VC_METRICS = [
    Metric(
        m.id,
        m.title,
        in_genstats=m.in_genstats,
        in_own_tabl=m.in_own_tabl,
        descr=m.descr,
        unit=m.unit,
        namespace=m.namespace or NAMESPACE,
        the_higher_the_worse=m.the_higher_the_worse,
    )
    for m in [
        # id_in_data                                        title (display name)   gen_stats  vc_table  unit description
        # Read stats:
        Metric("Total", "Variants", "#", "#", "", "Total number of variants (SNPs + MNPs + INDELS)."),
        Metric(
            "Biallelic",
            "Biallelic",
            None,
            "hid",
            "",
            "Number of sites in a genome that contains two observed alleles, counting the reference as one, and therefore allowing for one variant allele",
        ),
        Metric(
            "Multiallelic",
            "Multiallelic",
            "hid",
            "%",
            "",
            "Number of sites in the VCF that contain three or more observed alleles. The reference is counted as one, therefore allowing for two or more variant alleles",
        ),
        Metric(
            "SNPs",
            "SNP",
            "hid",
            "%",
            "",
            "Number of SNPs in the variant set. A variant is counted as an SNP when the reference, allele 1, and allele2 are all length 1",
        ),
        Metric("Indels", "Indel", "hid", "hid", "", "Number of insetions and deletions in the variant set."),
        Metric("Insertions", "Ins", None, "%", "", "Number of insetions in the variant set."),
        Metric("Deletions", "Del", None, "%", "", "Number of deletions in the variant set."),
        Metric(
            "Insertions (Hom)", "Hom ins", None, "hid", "", "Number of variants that contains homozygous insertions"
        ),
        Metric(
            "Insertions (Het)",
            "Het ins",
            None,
            "hid",
            "",
            "Number of variants where both alleles are insertions, but not homozygous",
        ),
        Metric("Deletions (Hom)", "Hom del", None, "hid", "", "Number of variants that contains homozygous deletions"),
        Metric(
            "Deletions (Het)",
            "Het del",
            None,
            "hid",
            "",
            "Number of variants where both alleles are deletion, but not homozygous",
        ),
        Metric(
            "Indels (Het)",
            "Het indel",
            None,
            "hid",
            "",
            "Number of variants where genotypes are either [insertion+deletion], [insertion+snp] or [deletion+snp].",
        ),
        Metric(
            "DeNovo SNPs",
            "DeNovo SNPs",
            None,
            None,
            "",
            "Number of DeNovo marked SNPs, with DQ > 0.05. Set the --qc-snp-denovo-quality-threshold option to the required threshold. The default is 0.05.",
        ),
        Metric(
            "DeNovo INDELs",
            "DeNovo indel",
            None,
            None,
            "",
            "Number of DeNovo marked indels, with DQ > 0.05. Set the --qc-snp-denovo-quality-threshold option to the required threshold. The default is 0.05.",
        ),
        Metric(
            "DeNovo MNPs",
            "DeNovo MNPs",
            None,
            None,
            "",
            "Number of DeNovo marked MNPs, with DQ > 0.05. Set the --qc-snp-denovo-quality-threshold option to the required threshold. The default is 0.05.",
        ),
        Metric(
            "Chr X number of SNPs over genome",
            "ChrX SNP",
            None,
            None,
            "",
            "Number of SNPs in chromosome X (or in the intersection of chromosome X with the target region). " "",
            "If there was no alignment to either chromosome X, this metric shows as NA",
        ),
        Metric(
            "Chr Y number of SNPs over genome",
            "ChrY SNP",
            None,
            None,
            "",
            "Number of SNPs in chromosome Y (or in the intersection of chromosome Y with the target region). " "",
            "If there was no alignment to either chromosome Y, this metric shows as NA",
        ),
        Metric(
            "(Chr X SNPs)/(chr Y SNPs) ratio over genome",
            "X/Y SNP ratio",
            None,
            "hid",
            "",
            "Number of SNPs in chromosome X (or in the intersection of chromosome X with the target region) " "",
            "divided by the number of SNPs in chromosome Y (or in the intersection of chromosome Y with the " "",
            "target region). If there was no alignment to either chromosome X or chromosome Y, this metric " "",
            "shows as NA",
        ),
        Metric(
            "SNP Transitions",
            "SNP Ti",
            None,
            None,
            "",
            "Number of transitions - interchanges of two purines (A<->G) or two pyrimidines (C<->T)",
        ),
        Metric(
            "SNP Transversions",
            "SNP Tv",
            None,
            None,
            "",
            "Number of transversions - interchanges of purine and pyrimidine bases",
        ),
        Metric("Ti/Tv ratio", "Ti/Tv", "hid", "#", "", "Ti/Tv ratio: ratio of transitions to transitions."),
        Metric("Heterozygous", "Het", "hid", "hid", "", "Number of heterozygous variants"),
        Metric("Homozygous", "Hom", "hid", "hid", "", "Number of homozygous variants"),
        Metric("Het/Hom ratio", "Het/Hom", "hid", "#", "", "Heterozygous/ homozygous ratio"),
        Metric(
            "In dbSNP",
            "In dbSNP",
            None,
            None,
            "",
            "Number of variants detected that are present in the dbsnp reference file. If no dbsnp file " "",
            "is provided via the --bsnp option, then both the In dbSNP and Novel metrics show as NA.",
            the_higher_the_worse=True,
        ),
        Metric(
            "Not in dbSNP",
            "Novel",
            None,
            None,
            "",
            "Number of all variants minus number of variants in dbSNP. If no dbsnp file " "",
            "is provided via the --bsnp option, then both the In dbSNP and Novel metrics show as NA.",
        ),
        Metric(
            "Percent Callability",
            "Callability",
            None,
            "#",
            "",
            "Available only in germline mode with gVCF output. The percentage of non-N reference " "",
            "positions having a PASSing genotype call. Multi-allelic variants are not counted. " "",
            "Deletions are counted for all the deleted reference positions only for homozygous calls. " "",
            "Only autosomes and chromosomes X, Y and M are considered.",
        ),
        Metric(
            "Percent Autosome Callability",
            "Autosome callability",
            None,
            "hid",
            "",
            "Available only in germline mode with gVCF output. The percentage of non-N reference " "",
            "positions having a PASSing genotype call. Multi-allelic variants are not counted. " "",
            "Deletions are counted for all the deleted reference positions only for homozygous calls. " "",
            "Only autosomes are considered (for all chromosomes, see the Callability metric).",
        ),
        Metric(
            "Filtered vars",
            "Filt var",
            "hid",
            "hid",
            "",
            "Number of raw variants minus the number of PASSed variants",
            the_higher_the_worse=True,
        ),
        Metric(
            "Filtered SNPs",
            "Filt SNP",
            "hid",
            "%",
            "",
            "Number of raw SNPs minus the number of PASSed SNPs",
            the_higher_the_worse=True,
        ),
        Metric(
            "Filtered indels",
            "Filt indel",
            "hid",
            "%",
            "",
            "Number of raw indels minus the number of PASSed indels",
            the_higher_the_worse=True,
        ),
        Metric(
            "Reads Processed",
            "VC reads",
            None,
            "#",
            "reads",
            "The number of reads used for variant calling, excluding any duplicate marked reads and reads falling outside of the target region",
        ),
    ]
]


def parse_vc_metrics_file(f):
    """
    T_SRR7890936_50pc.vc_metrics.csv or T_SRR7890936_50pc.gvcf_metrics.csv

    VARIANT CALLER SUMMARY,,Number of samples,1
    VARIANT CALLER SUMMARY,,Reads Processed,2721782043
    VARIANT CALLER SUMMARY,,Child Sample,NA
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Total,170013,100.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Biallelic,170013,100.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Multiallelic,0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,SNPs,138978,81.75
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Insertions (Hom),0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Insertions (Het),14291,8.41
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Deletions (Hom),0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Deletions (Het),16744,9.85
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Indels (Het),0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Chr X number of SNPs over genome,32487
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Chr Y number of SNPs over genome,131
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,(Chr X SNPs)/(chr Y SNPs) ratio over genome,247.99
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,SNP Transitions,79948
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,SNP Transversions,59025
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Ti/Tv ratio,1.35
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Heterozygous,170013
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Homozygous,0
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Het/Hom ratio,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,In dbSNP,0,0.00
    VARIANT CALLER PREFILTER,T_SRR7890936_50pc,Not in dbSNP,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Total,123219,100.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Biallelic,123219,100.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Multiallelic,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,SNPs,104900,85.13
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Insertions (Hom),0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Insertions (Het),8060,6.54
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Deletions (Hom),0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Deletions (Het),10259,8.33
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Indels (Het),0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Chr X number of SNPs over genome,28162
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Chr Y number of SNPs over genome,8
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,(Chr X SNPs)/(chr Y SNPs) ratio over genome,3520.25
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,SNP Transitions,62111
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,SNP Transversions,42789
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Ti/Tv ratio,1.45
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Heterozygous,123219
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Homozygous,0
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Het/Hom ratio,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,In dbSNP,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Not in dbSNP,0,0.00
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Percent Callability,NA
    VARIANT CALLER POSTFILTER,T_SRR7890936_50pc,Percent Autosome Callability,NA
    """

    summary_data = dict()
    prefilter_data = dict()
    postfilter_data = dict()

    for line in f["f"].splitlines():
        fields = line.split(",")
        analysis = fields[0]
        # sample = fields[1]
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

        if analysis == "VARIANT CALLER SUMMARY":
            summary_data[metric] = value

        if analysis == "VARIANT CALLER PREFILTER":
            prefilter_data[metric] = value

        if analysis in ["VARIANT CALLER POSTFILTER", "VARIANT CALLER POSTFILTER GVCF"]:
            postfilter_data[metric] = value
            if percentage is not None:
                postfilter_data[metric + " pct"] = percentage

    # adding few more metrics: total insertions, deletions and indels numbers
    for data in [prefilter_data, postfilter_data]:
        if exist_and_number(data, "Insertions (Hom)", "Insertions (Het)"):
            data["Insertions"] = data["Insertions (Hom)"] + data["Insertions (Het)"]

        if exist_and_number(data, "Deletions (Hom)", "Deletions (Het)"):
            data["Deletions"] = data["Deletions (Hom)"] + data["Deletions (Het)"]

        if exist_and_number(data, "Insertions", "Deletions"):
            data["Indels"] = data["Insertions"] + data["Deletions"]

        if exist_and_number(data, "Total") and data["Total"] != 0:
            if exist_and_number(data, "Insertions"):
                data["Insertions pct"] = data["Insertions"] / data["Total"] * 100.0
            if exist_and_number(data, "Deletions"):
                data["Deletions pct"] = data["Deletions"] / data["Total"] * 100.0
            if exist_and_number(data, "Indels"):
                data["Indels pct"] = data["Indels"] / data["Total"] * 100.0

    data = postfilter_data
    data.update(summary_data)

    # we are not really interested in all the details of pre-filtered variants, however
    # it would be nice to report how much we filtered out
    if exist_and_number(data, "Total") and exist_and_number(prefilter_data, "Total"):
        # prefilter numbers can be zero (see https://github.com/MultiQC/MultiQC/issues/2611), so only reporting
        # if they are non-zero, otherwise will get negative number for "Filtered vars"
        if prefilter_data["Total"] > 0:
            data["Filtered vars"] = prefilter_data["Total"] - data["Total"]

    if exist_and_number(data, "SNPs") and exist_and_number(prefilter_data, "SNPs"):
        if prefilter_data["SNPs"] > 0:
            data["Filtered SNPs"] = prefilter_data["SNPs"] - data["SNPs"]

    if exist_and_number(data, "Indels") and exist_and_number(prefilter_data, "Indels"):
        if prefilter_data["Indels"] > 0:
            data["Filtered indels"] = prefilter_data["Indels"] - data["Indels"]

    if (
        exist_and_number(prefilter_data, "Total")
        and exist_and_number(data, "Filtered vars")
        and prefilter_data["Total"] != 0
    ):
        data["Filtered vars pct"] = data["Filtered vars"] / prefilter_data["Total"] * 100.0

    if (
        exist_and_number(prefilter_data, "SNPs")
        and exist_and_number(data, "Filtered SNPs")
        and prefilter_data["SNPs"] != 0
    ):
        data["Filtered SNPs pct"] = data["Filtered SNPs"] / prefilter_data["SNPs"] * 100.0

    if (
        exist_and_number(prefilter_data, "Indels")
        and exist_and_number(data, "Filtered Indels")
        and prefilter_data["Indels"] != 0
    ):
        data["Filtered indels pct"] = data["Filtered indels"] / prefilter_data["Indels"] * 100.0

    return data
