# coding=utf-8


import itertools
import logging
from collections import defaultdict

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

from .utils import Metric, exist_and_number, make_headers

# Initialise the logger
log = logging.getLogger(__name__)


NAMESPACE = "Mapping"


class DragenMappingMetics(BaseMultiqcModule):
    def add_mapping_metrics(self):
        data_by_rg_by_sample = defaultdict(dict)
        data_by_phenotype_by_sample = defaultdict(dict)

        for f in self.find_log_files("dragen/mapping_metrics"):
            data_by_readgroup, data_by_phenotype = parse_mapping_metrics_file(f)
            s_name = f["s_name"]
            if s_name in data_by_rg_by_sample:
                log.debug(f"Duplicate DRAGEN output prefix found! Overwriting: {s_name}")
            self.add_data_source(f, section="mapping_metrics")
            data_by_phenotype_by_sample[s_name].update(data_by_phenotype)

            for phenotype, phenotype_d in data_by_readgroup.items():
                for rg, data in phenotype_d.items():
                    if any(rg in d_rg for sn, d_rg in data_by_rg_by_sample.items()):
                        log.debug(f"Duplicate read group name {rg} found for output prefix {s_name}! Overwriting")
            data_by_rg_by_sample[s_name].update(data_by_readgroup)

            # Superfluous function call to confirm that it is used in this module
            # Replace None with actual version if it is available
            self.add_software_version(None, s_name)

        # filter to strip out ignored sample names:
        data_by_rg_by_sample = self.ignore_samples(data_by_rg_by_sample)
        data_by_phenotype_by_sample = self.ignore_samples(data_by_phenotype_by_sample)

        # flattening phenotype-sample data by adding a prefix " normal" to the normal samples
        data_by_sample = dict()
        for sn in data_by_phenotype_by_sample:
            for phenotype in data_by_phenotype_by_sample[sn]:
                new_sn = sn
                if phenotype == "normal":
                    new_sn = sn + "_normal"
                data_by_sample[new_sn] = data_by_phenotype_by_sample[sn][phenotype]

        # flattening phenotype-sample data by adding a prefix " normal" to the normal samples
        data_by_rg_by_sample_new = dict()
        for sn in data_by_rg_by_sample:
            for phenotype in data_by_rg_by_sample[sn]:
                new_sn = sn
                if phenotype == "normal":
                    new_sn = sn + "_normal"
                data_by_rg_by_sample_new[new_sn] = data_by_rg_by_sample[sn][phenotype]
        data_by_rg_by_sample = data_by_rg_by_sample_new
        del data_by_rg_by_sample_new

        if not data_by_rg_by_sample and not data_by_phenotype_by_sample:
            return set()

        # Write data to file
        self.write_data_file(data_by_sample, "dragen_map_metrics")

        self.report_mapping_metrics(data_by_sample, data_by_rg_by_sample)

        return data_by_rg_by_sample.keys()

    def report_mapping_metrics(self, data_by_sample, data_by_rg_by_sample):
        # merging all read group data
        data_by_rg = dict()

        for sname in data_by_rg_by_sample:
            for rg, d in data_by_rg_by_sample[sname].items():
                # Check only one read-group per sample
                if len(list(data_by_rg_by_sample[sname].keys())) == 1:
                    rg = sname
                    data_by_rg[rg] = d
                    continue

                # Multiple read-groups per sample
                rg = sname + f" ({rg}) "
                data_by_rg[rg] = d

        # getting all available metric names to determine table headers
        all_metric_names = set()
        for sn, d_by_rg in data_by_rg_by_sample.items():
            for rg, data in d_by_rg.items():
                for m in data.keys():
                    all_metric_names.add(m)
        # and making headers
        genstats_headers, own_tabl_headers = make_headers(all_metric_names, MAPPING_METRICS)

        self.general_stats_addcols(data_by_sample, genstats_headers, namespace=NAMESPACE)

        self.add_section(
            name="Mapping metrics",
            anchor="dragen-mapping-metrics",
            description="""
            Mapping metrics, similar to the metrics computed by the samtools-stats command.
            Shown on per read group level. To see per-sample level metrics, refer to the general
            stats table.
            """,
            plot=table.plot(
                data_by_rg,
                own_tabl_headers,
                pconfig={
                    "id": "dragen-mapping-metrics-table",
                    "namespace": NAMESPACE,
                    "title": "DRAGEN: Mapping metrics",
                },
            ),
        )

        # Skip adding the barplot if it's not informative, such as if all
        #  reads are unmapped due to a FastQcOnly workflow
        all_unmapped = True
        unmapped_key = "Unmapped reads pct"
        for rg, data in data_by_rg.items():
            if unmapped_key in data and data[unmapped_key] < 100.0:
                all_unmapped = False

        if all_unmapped:
            return set()

        # Make bargraph plots of mapped, dupped and paired reads
        self.__map_pair_dup_read_chart(data_by_rg)

    def __map_pair_dup_read_chart(self, data_by_sample):
        paired_reads_data = {}
        mapped_reads_data = {}
        category_labels = []
        data_labels = []

        add_paired_label = False
        add_mapped_label = False
        for sample_id, data in data_by_sample.items():
            # Dragen 3.9 has replaced 'rRNA filtered reads' with 'Adjustment of reads matching filter contigs'
            if "rRNA filtered reads" in data.keys():
                rrna_filtered_reads_key = "rRNA filtered reads"
            elif "Adjustment of reads matching filter contigs" in data.keys():
                rrna_filtered_reads_key = "Adjustment of reads matching filter contigs"
            else:
                rrna_filtered_reads_key = None

            if data.get("Mapped reads R2", None) == 0:
                log.debug(f"single-ended data detected, skipping mapping/paired percentages plot for: {sample_id}")
            elif data.get("Not properly paired reads (discordant)", 0) + data.get(
                "Properly paired reads", 0
            ) + data.get("Singleton reads (itself mapped; mate unmapped)", 0) + data.get("Unmapped reads", 0) + (
                data.get(rrna_filtered_reads_key, 0)
                if rrna_filtered_reads_key is not None and rrna_filtered_reads_key == "rRNA filtered reads"
                else 0
            ) != data.get("Total reads in RG", 0):
                log.warning(
                    "sum of unpaired/discordant/proppaired/unmapped reads not matching total, "
                    "skipping mapping/paired percentages plot for: {}".format(sample_id)
                )
            else:
                paired_reads_data[sample_id] = data
                add_paired_label = True

            if data.get("Number of unique & mapped reads (excl. duplicate marked reads)", 0) + data.get(
                "Number of duplicate marked reads", 0
            ) + data.get("Unmapped reads", 0) + (
                data.get(rrna_filtered_reads_key, 0)
                if rrna_filtered_reads_key is not None and rrna_filtered_reads_key == "rRNA filtered reads"
                else 0
            ) != data.get("Total reads in RG", 0):
                log.warning(
                    "sum of unique/duplicate/unmapped reads not matching total, "
                    "skipping mapping/duplicates percentages plot for: {}".format(sample_id)
                )
            else:
                mapped_reads_data[sample_id] = data
                add_mapped_label = True

        # Add labels
        if add_paired_label:
            category_labels.append(
                {
                    "Properly paired reads": {"color": "#099109", "name": "Paired, properly"},
                    "Not properly paired reads (discordant)": {"color": "#c27a0e", "name": "Paired, discordant"},
                    "Singleton reads (itself mapped; mate unmapped)": {"color": "#912476", "name": "Singleton"},
                    "Unmapped reads": {"color": "#b1084c", "name": "Unmapped"},
                }
            )
            data_labels.append(
                {"name": "Paired vs. discordant vs. singleton", "ylab": "Reads", "cpswitch_counts_label": "Reads"}
            )
        if add_mapped_label:
            mapped_chart_labels = {
                "Number of unique & mapped reads (excl. duplicate marked reads)": {
                    "color": "#437bb1",
                    "name": "Unique",
                },
                "Number of duplicate marked reads": {"color": "#f5a742", "name": "Duplicated"},
                "Unmapped reads": {"color": "#b1084c", "name": "Unmapped"},
            }
            if "Adjustment of reads matching filter contigs" in next(iter(data_by_sample.values())):
                mapped_chart_labels["Adjustment of reads matching filter contigs"] = {
                    "color": "#43b14a",
                    "name": "rRNA filtered",
                }
            elif "rRNA filtered reads" in next(iter(data_by_sample.values())):
                mapped_chart_labels["rRNA filtered reads"] = {"color": "#43b14a", "name": "rRNA filtered"}
            elif "Adjustment of reads matching filter contigs" in next(iter(data_by_sample.values())):
                mapped_chart_labels["Adjustment of reads matching filter contigs"] = {
                    "color": "#43b14a",
                    "name": "rRNA filtered",
                }
            category_labels.append(mapped_chart_labels)
            data_labels.append(
                {"name": "Unique vs duplicated vs unmapped", "ylab": "Reads", "cpswitch_counts_label": "Reads"}
            )

        data_to_plot = [d for d in [paired_reads_data, mapped_reads_data] if d]  # Leaves out empty dicts
        self.add_section(
            name="Mapped / paired / duplicated",
            anchor="dragen-mapped-paired-duplicated",
            description="Distribution of reads based on pairing, duplication and mapping.",
            plot=bargraph.plot(
                data_to_plot,
                category_labels,
                {
                    "id": "mapping_dup_percentage_plot",
                    "title": "Dragen: Mapped/paired/duplicated reads per read group",
                    "ylab": "Reads",
                    "cpswitch_counts_label": "Reads",
                    "data_labels": data_labels,
                },
            ),
        )


def parse_mapping_metrics_file(f):
    """
    Mapping and aligning metrics, like the metrics computed by the Samtools Flagstat command, are available
    on an aggregate level (over all input data), and on a per read group level. Unless explicitly stated,
    the metrics units are in reads (ie, not in terms of pairs or alignments).

    T_SRR7890936_50pc.mapping_metrics.csv

    # phenotype-level metrics (tumor or normal):
    TUMOR MAPPING/ALIGNING SUMMARY,,Total input reads,2200000000,100.00
    TUMOR MAPPING/ALIGNING SUMMARY,,Number of duplicate marked reads,433637413,19.71
    TUMOR MAPPING/ALIGNING SUMMARY,,Number of duplicate marked and mate reads removed,NA
    TUMOR MAPPING/ALIGNING SUMMARY,,Number of unique reads (excl. duplicate marked reads),1766362587,80.29
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with mate sequenced,2200000000,100.00
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads without mate sequenced,0,0.00
    TUMOR MAPPING/ALIGNING SUMMARY,,QC-failed reads,0,0.00
    TUMOR MAPPING/ALIGNING SUMMARY,,Mapped reads,2130883930,96.86
    TUMOR MAPPING/ALIGNING SUMMARY,,Mapped reads R1,1066701794,96.97
    TUMOR MAPPING/ALIGNING SUMMARY,,Mapped reads R2,1064182136,96.74
    TUMOR MAPPING/ALIGNING SUMMARY,,Number of unique & mapped reads (excl. duplicate marked reads),1697246517,77.15
    TUMOR MAPPING/ALIGNING SUMMARY,,Unmapped reads,69116070,3.14
    TUMOR MAPPING/ALIGNING SUMMARY,,Singleton reads (itself mapped; mate unmapped),3917092,0.18
    TUMOR MAPPING/ALIGNING SUMMARY,,Paired reads (itself & mate mapped),2126966838,96.68
    TUMOR MAPPING/ALIGNING SUMMARY,,Properly paired reads,2103060370,95.59
    TUMOR MAPPING/ALIGNING SUMMARY,,Not properly paired reads (discordant),23906468,1.09
    TUMOR MAPPING/ALIGNING SUMMARY,,Paired reads mapped to different chromosomes,17454370,0.82
    TUMOR MAPPING/ALIGNING SUMMARY,,Paired reads mapped to different chromosomes (MAPQ>=10),6463547,0.30
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [40:inf),2002661377,91.03
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [30:40),7169392,0.33
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [20:30),16644390,0.76
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [10:20),20280057,0.92
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with MAPQ [ 0:10),84128714,3.82
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with MAPQ NA (Unmapped reads),69116070,3.14
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with indel R1,26849051,2.52
    TUMOR MAPPING/ALIGNING SUMMARY,,Reads with indel R2,24810803,2.33
    TUMOR MAPPING/ALIGNING SUMMARY,,Total bases,330000000000
    TUMOR MAPPING/ALIGNING SUMMARY,,Total bases R1,165000000000
    TUMOR MAPPING/ALIGNING SUMMARY,,Total bases R2,165000000000
    TUMOR MAPPING/ALIGNING SUMMARY,,Mapped bases R1,160005269100
    TUMOR MAPPING/ALIGNING SUMMARY,,Mapped bases R2,159627320400
    TUMOR MAPPING/ALIGNING SUMMARY,,Soft-clipped bases R1,1757128997,1.10
    TUMOR MAPPING/ALIGNING SUMMARY,,Soft-clipped bases R2,3208748350,2.01
    TUMOR MAPPING/ALIGNING SUMMARY,,Mismatched bases R1,585802788,0.37
    TUMOR MAPPING/ALIGNING SUMMARY,,Mismatched bases R2,1155805091,0.72
    TUMOR MAPPING/ALIGNING SUMMARY,,Mismatched bases R1 (excl. indels),501394281,0.31
    TUMOR MAPPING/ALIGNING SUMMARY,,Mismatched bases R2 (excl. indels),1073788605,0.67
    TUMOR MAPPING/ALIGNING SUMMARY,,Q30 bases,297564555927,90.17
    TUMOR MAPPING/ALIGNING SUMMARY,,Q30 bases R1,155492239719,94.24
    TUMOR MAPPING/ALIGNING SUMMARY,,Q30 bases R2,142072316208,86.10
    TUMOR MAPPING/ALIGNING SUMMARY,,Q30 bases (excl. dups & clipped bases),246555769158
    TUMOR MAPPING/ALIGNING SUMMARY,,Total alignments,2190085267
    TUMOR MAPPING/ALIGNING SUMMARY,,Secondary alignments,0
    TUMOR MAPPING/ALIGNING SUMMARY,,Supplementary (chimeric) alignments,59201337
    TUMOR MAPPING/ALIGNING SUMMARY,,Estimated read length,150.00
    TUMOR MAPPING/ALIGNING SUMMARY,,Average sequenced coverage over genome,102.83
    TUMOR MAPPING/ALIGNING SUMMARY,,Insert length: mean,383.00
    TUMOR MAPPING/ALIGNING SUMMARY,,Insert length: median,376.00
    TUMOR MAPPING/ALIGNING SUMMARY,,Insert length: standard deviation,85.15
    # general metrics - reporting in the header (except for DRAGEN mapping rate may be different to T and N:
    TUMOR MAPPING/ALIGNING SUMMARY,,Bases in reference genome,3209286105
    TUMOR MAPPING/ALIGNING SUMMARY,,Bases in target bed [% of genome],NA
    TUMOR MAPPING/ALIGNING SUMMARY,,Provided sex chromosome ploidy,NA
    TUMOR MAPPING/ALIGNING SUMMARY,,DRAGEN mapping rate [mil. reads/second],0.39
    # then same for normal:
    NORMAL MAPPING/ALIGNING SUMMARY,,Total input reads,1100000000,100.00
    NORMAL MAPPING/ALIGNING SUMMARY,,Number of duplicate marked reads,123518125,11.23
    ...
    # then tumor and normal per-read-group metrics::
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Total reads in RG,2200000000,100.00
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Number of duplicate marked reads,433637413,19.71
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Number of duplicate marked and mate reads removed,NA
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Number of unique reads (excl. duplicate marked reads),1766362587,80.29
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with mate sequenced,2200000000,100.00
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads without mate sequenced,0,0.00
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,QC-failed reads,0,0.00
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mapped reads,2130883930,96.86
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mapped reads R1,1066701794,96.97
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mapped reads R2,1064182136,96.74
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Number of unique & mapped reads (excl. duplicate marked reads),1697246517,77.15
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Unmapped reads,69116070,3.14
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Singleton reads (itself mapped; mate unmapped),3917092,0.18
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Paired reads (itself & mate mapped),2126966838,96.68
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Properly paired reads,2103060370,95.59
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Not properly paired reads (discordant),23906468,1.09
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Paired reads mapped to different chromosomes,17454370,0.82
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Paired reads mapped to different chromosomes (MAPQ>=10),6463547,0.30
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with MAPQ [40:inf),2002661377,91.03
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with MAPQ [30:40),7169392,0.33
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with MAPQ [20:30),16644390,0.76
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with MAPQ [10:20),20280057,0.92
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with MAPQ [ 0:10),84128714,3.82
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with MAPQ NA (Unmapped reads),69116070,3.14
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with indel R1,26849051,2.52
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Reads with indel R2,24810803,2.33
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Total bases,330000000000
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Total bases R1,165000000000
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Total bases R2,165000000000
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mapped bases R1,160005269100
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mapped bases R2,159627320400
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Soft-clipped bases R1,1757128997,1.10
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Soft-clipped bases R2,3208748350,2.01
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mismatched bases R1,585802788,0.37
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mismatched bases R2,1155805091,0.72
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mismatched bases R1 (excl. indels),501394281,0.31
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Mismatched bases R2 (excl. indels),1073788605,0.67
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Q30 bases,297564555927,90.17
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Q30 bases R1,155492239719,94.24
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Q30 bases R2,142072316208,86.10
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Q30 bases (excl. dups & clipped bases),246555769158
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Total alignments,2190085267
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Secondary alignments,0
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Supplementary (chimeric) alignments,59201337
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Estimated read length,150.00
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Average sequenced coverage over genome,102.83
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Insert length: mean,383.01
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Insert length: median,376.00
    TUMOR MAPPING/ALIGNING PER RG,T_SRR7890936_50pc,Insert length: standard deviation,87.58
    # same for normal:
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Total reads in RG,1100000000,100.00
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of duplicate marked reads,NA
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of duplicate marked and mate reads removed,NA
    NORMAL MAPPING/ALIGNING PER RG,N_SRR7890889,Number of unique reads (excl. duplicate marked reads),NA
    ...

    We are reporting summary metrics in the general stats table, and per-read-group in a separate table.
    """

    data_by_readgroup = defaultdict(dict)
    data_by_phenotype = defaultdict(dict)

    for line in f["f"].splitlines():
        fields = line.split(",")
        phenotype = fields[0].split("/")[0].split(" ")[0].lower()  # TUMOR MAPPING -> tumor
        analysis = fields[0].split("/")[1]  # ALIGNING SUMMARY, ALIGNING PER RG
        metric = fields[2]
        value = fields[3] if fields[3] != "NA" else None
        try:
            value = int(value)
        except (ValueError, TypeError):
            try:
                value = float(value)
            except (ValueError, TypeError):
                pass

        percentage = None
        if len(fields) > 4:  # percentage
            percentage = fields[4]
            try:
                percentage = float(percentage)
            except ValueError:
                pass

        # sample-unspecific metrics are reported only in ALIGNING SUMMARY sections
        if analysis == "ALIGNING SUMMARY":
            if value is not None:
                data_by_phenotype[phenotype][metric] = value
            if percentage is not None:
                data_by_phenotype[phenotype][metric + " pct"] = percentage

        # for sample-specific metrics, using ALIGNING PER RG because it has the sample name in the 2nd col
        if analysis == "ALIGNING PER RG":
            # setting normal and tumor sample names for future use
            readgroup = fields[1]
            if readgroup not in data_by_readgroup[phenotype].keys():
                data_by_readgroup[phenotype][readgroup] = {}
            if value is not None:
                data_by_readgroup[phenotype][readgroup][metric] = value
            if percentage is not None:
                data_by_readgroup[phenotype][readgroup][metric + " pct"] = percentage

    # adding some missing values that we wanna report for consistency
    # Expand data_by_readgroup to values below phenotype level
    for data in itertools.chain(
        *[data_by_readgroup[key].values() for key in data_by_readgroup.keys()], data_by_phenotype.values()
    ):
        # fixing when deduplication wasn't performed, or running with single-end data
        for field in [
            "Number of duplicate marked and mate reads removed",
            "Number of unique reads (excl. duplicate marked reads)",
            "Mismatched bases R2 (excl. indels)",
            "Mismatched bases R2",
            "Soft-clipped bases R2",
        ]:
            try:
                data.pop(field)
            except KeyError:
                pass

        # adding alignment percentages
        if exist_and_number(data, "Total alignments", "Secondary alignments") and data["Total alignments"] > 0:
            data["Secondary alignments pct"] = data["Secondary alignments"] / data["Total alignments"] * 100.0

        # adding some missing bases percentages
        if exist_and_number(data, "Total bases") and data["Total bases"] > 0:
            for m in ["Q30 bases (excl. dups & clipped bases)", "Mapped bases R1", "Mapped bases R2"]:
                if exist_and_number(data, m):
                    data[m + " pct"] = data[m] / data["Total bases"] * 100.0

    return data_by_readgroup, data_by_phenotype


MAPPING_METRICS = [
    # id_in_data                                                    # title (display name)  # in_genstats  # in_own_tabl # unit  # description
    # Read stats:
    Metric(
        "Total input reads",
        "Input reads",
        "#",
        None,
        "reads",
        "Total number of input reads for this sample (or total number of reads in all input read groups combined), {}",
    ),
    Metric("Total reads in RG", "Input reads", None, "#", "reads", "Total number of reads in this RG, {}"),
    Metric(
        "Average sequenced coverage over genome",
        "Raw depth",
        "hid",
        "#",
        "x",
        "Average sequenced coverage over genome, based on all input reads (including duplicate, clipped and low quality bases and reads)",
    ),
    Metric("Reads with mate sequenced", "Paired", "hid", "%", "reads", "Number of reads with a mate sequenced, {}"),
    Metric(
        "Reads without mate sequenced",
        "Single",
        None,
        "hid",
        "reads",
        "Number of reads without a mate sequenced, {}",
        the_higher_the_worse=True,
    ),
    Metric(
        "QC-failed reads",
        "QC-fail",
        "hid",
        "%",
        "reads",
        "Number of reads not passing platform/vendor quality checks (SAM flag 0x200), {}",
        precision=2,
        the_higher_the_worse=True,
    ),
    Metric("Mapped reads", "Map", "hid", "hid", "reads", "Number of mapped reads, {}"),
    Metric("Mapped reads R1", "Map R1", None, "hid", "reads", "Number of mapped reads R1, {}"),
    Metric("Mapped reads R2", "Map R2", None, "hid", "reads", "Number of mapped reads R2, {}"),
    Metric("Unmapped reads", "Unmap", "%", "%", "reads", "Number of unmapped reads, {}", the_higher_the_worse=True),
    Metric(
        "rRNA filtered reads", "rRNA", "%", "%", "reads", "Number of rRNA filtered reads, {}", the_higher_the_worse=True
    ),
    Metric("Reads with MAPQ [40:inf)", "MQ⩾40", None, "hid", "reads", "Number of reads with MAPQ [40:inf), {}"),
    Metric(
        "Number of duplicate marked reads",
        "Dup",
        "%",
        "%",
        "reads",
        "Number of duplicate marked reads as a result of the --enable-duplicatemarking option being used, {}",
        the_higher_the_worse=True,
    ),
    Metric(
        "Number of duplicate marked and mate reads removed",
        "Dup with mates removed",
        None,
        "hid",
        "reads",
        "Number of reads marked as duplicates, along with any mate reads, that are removed "
        "when the --remove-duplicates option is used, {}",
        the_higher_the_worse=True,
    ),
    Metric(
        "Number of unique reads (excl. duplicate marked reads)",
        "Uniq",
        None,
        "hid",
        "reads",
        "Number of unique reads (all reads minus duplicate marked reads), {}",
    ),
    Metric(
        "Number of unique & mapped reads (excl. duplicate marked reads)",
        "Uniq map",
        "hid",
        "hid",
        "reads",
        "Number of unique & mapped reads (mapped reads minus duplicate marked reads), {}",
    ),
    Metric(
        "Paired reads (itself & mate mapped)",
        "Self & mate map",
        None,
        "hid",
        "reads",
        "Number of reads mapped in pairs (itself & mate mapped), {}",
    ),
    Metric(
        "Properly paired reads",
        "Prop pair",
        "%",
        "%",
        "reads",
        "Number of properly paired reads, {} (both reads in pair are mapped and "
        "fall within an acceptable range from each other based on the estimated insert length distribution). "
        "Duplicate and low quality reads are not excluded, i.e. the percentage is calculated relative to the total number of reads.",
    ),
    Metric(
        "Not properly paired reads (discordant)",
        "Discord",
        "hid",
        "%",
        "reads",
        "Number of discordant reads: paired reads minus properly paired reads , {}",
        precision=2,
        the_higher_the_worse=True,
    ),
    Metric(
        "Singleton reads (itself mapped; mate unmapped)",
        "Singleton",
        None,
        "%",
        "reads",
        "Number of singleton reads: itself mapped; mate unmapped, {}",
        precision=2,
        the_higher_the_worse=True,
    ),
    Metric(
        "Paired reads mapped to different chromosomes",
        "Diff chrom",
        None,
        "hid",
        "reads",
        "Number of paired reads with a mate mapped to a different chromosome, {}",
        precision=2,
        the_higher_the_worse=True,
    ),
    Metric(
        "Paired reads mapped to different chromosomes (MAPQ>=10)",
        "Diff chr, MQ⩾10",
        "hid",
        "%",
        "reads",
        "Number of paired reads, mapped with MAPQ>=10 and with a mate mapped to a different chromosome, {}",
        precision=2,
        the_higher_the_worse=True,
    ),
    Metric(
        "Reads with indel R1",
        "Reads with indel R1",
        None,
        "hid",
        "reads",
        "Number of R1 reads containing at least 1 indel, {}",
        the_higher_the_worse=True,
    ),
    Metric(
        "Reads with indel R2",
        "Reads with indel R2",
        None,
        "hid",
        "reads",
        "Number of R2 reads containing at least 1 indel, {}",
        the_higher_the_worse=True,
    ),
    # Read length stats
    Metric(
        "Estimated read length",
        "Read len",
        "hid",
        "hid",
        "bp",
        "Estimated read length. Total number of input bases divided by the number of reads",
    ),
    Metric("Insert length: mean", "Avg IS", "hid", "hid", "bp", "Mean insert size"),
    Metric("Insert length: median", "Med IS", "#", "#", "bp", "Median insert size", precision=0),
    Metric(
        "Insert length: standard deviation", "IS std", "hid", "hid", "bp", "Standard deviation of insert size deviation"
    ),
    # Bases stats:
    Metric("Total bases", "Input bases", "hid", "hid", "bases", "Total number of bases sequenced, {}"),
    Metric("Total bases R1", "Input bases R1", None, "hid", "bases", "Total number of bases sequenced on R1 reads, {}"),
    Metric("Total bases R2", "Input bases R2", None, "hid", "bases", "Total number of bases sequenced on R2 reads, {}"),
    Metric("Mapped bases R1", "Mapped bases R1", None, "hid", "bases", "Number of mapped bases on R1 reads, {}"),
    Metric("Mapped bases R2", "Mapped bases R2", None, "hid", "bases", "Number of mapped bases on R2 reads, {}"),
    Metric(
        "Soft-clipped bases R1",
        "Soft-clip bases R1",
        None,
        "hid",
        "bases",
        "Number of soft-clipped bases on R1 reads, {}",
        the_higher_the_worse=True,
    ),
    Metric(
        "Soft-clipped bases R2",
        "Soft-clip bases R2",
        None,
        "hid",
        "bases",
        "Number of soft-clipped bases on R2 reads, {}",
        the_higher_the_worse=True,
    ),
    Metric(
        "Mismatched bases R1",
        "MM bases R1",
        None,
        "hid",
        "bases",
        "Number of mismatched bases on R1, {}, which is the sum of SNP count and indel lengths. "
        "It does not count anything within soft clipping, or RNA introns. "
        "It also does not count a mismatch if either the reference base or read base is N",
        the_higher_the_worse=True,
    ),
    Metric(
        "Mismatched bases R2",
        "MM bases R2",
        None,
        "hid",
        "bases",
        "Number of mismatched bases on R2, {}, which is the sum of SNP count and indel lengths. "
        "It does not count anything within soft clipping, or RNA introns. "
        "It also does not count a mismatch if either the reference base or read base is N",
        the_higher_the_worse=True,
    ),
    Metric(
        "Mismatched bases R1 (excl. indels)",
        "MM bases R1 excl indels",
        None,
        "hid",
        "bases",
        "Number of mismatched bases on R1, {}. The indels lengts are ignored. "
        "It does not count anything within soft clipping, or RNA introns. "
        "It also does not count a mismatch if either the reference base or read base is N",
        the_higher_the_worse=True,
    ),
    Metric(
        "Mismatched bases R2 (excl. indels)",
        "MM bases R2 excl indels",
        None,
        "hid",
        "bases",
        "Number of mismatched bases on R2, {}. The indels lengts are ignored. "
        "It does not count anything within soft clipping, or RNA introns. "
        "It also does not count a mismatch if either the reference base or read base is N",
        the_higher_the_worse=True,
    ),
    Metric("Q30 bases", "Q30", "hid", "hid", "bases", "Number of raw bases with BQ >= 30, {}"),
    Metric("Q30 bases R1", "Q30 R1", None, "hid", "bases", "Number of raw bases on R1 reads with BQ >= 30, {}"),
    Metric("Q30 bases R2", "Q30 R2", None, "hid", "bases", "Number of raw bases on R2 reads with BQ >= 30, {}"),
    Metric(
        "Q30 bases (excl. dups & clipped bases)",
        "Q30 excl dup & clipped",
        None,
        "hid",
        "bases",
        "Number of non-clipped bases with BQ >= 30 on non-duplicate reads, {}",
    ),
    # Alignments stats:
    Metric("Total alignments", "Alignments", "hid", "#", "reads", "Total number of alignments with MQ > 0, {}"),
    Metric(
        "Secondary alignments",
        "Sec'ry",
        "hid",
        "%",
        "reads",
        "Number of secondary alignments, {}. Secondary alignment occurs when "
        "a given read could align reasonably well to more than one place. "
        'One of the possible reported alignments is termed "primary" and '
        'the others will be marked as "secondary".',
        precision=2,
    ),
    Metric(
        "Supplementary (chimeric) alignments",
        "Suppl'ry",
        "hid",
        "%",
        "reads",
        "Number of supplementary (chimeric) alignments, {}. A chimeric read is split "
        "over multiple loci (possibly due to structural variants). One alignment is "
        "referred to as the representative alignment, the other are supplementary",
        precision=2,
    ),
    Metric(
        "Estimated sample contamination",
        "Contam'n",
        "#",
        "#",
        "proportion",
        "Fraction of estimated sample contamination",
        precision=2,
        the_higher_the_worse=True,
    ),
    Metric(
        "rRNA filtered reads", "rRNA", "%", "%", "reads", "Number of rRNA filtered reads, {}", the_higher_the_worse=True
    ),
    Metric(
        "Adjustment of reads matching filter contigs",
        "rRNA / Filtered Contigs",
        "%",
        "%",
        "reads",
        "Number of filtered reads, {}",
        the_higher_the_worse=True,
    ),
]
