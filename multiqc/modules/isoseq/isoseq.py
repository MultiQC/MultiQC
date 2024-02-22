from collections import defaultdict

import csv

import json

import logging

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, box

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Iso-Seq",
            anchor="isoseq",
            href="https://github.com/PacificBiosciences/IsoSeq",
            info="contains the newest tools to identify transcripts in PacBio single-molecule sequencing data (HiFi reads).",
            # doi=,  # Not published
        )

        refine_json_data_by_sample, refine_csv_data_by_sample = self._parse_refine()
        cnt_by_cluster_id_by_sample = self._parse_cluster()

        if not cnt_by_cluster_id_by_sample and not refine_json_data_by_sample and not refine_csv_data_by_sample:
            raise UserWarning

        # Superfluous function call to confirm that it is used in this module
        # Replace None with actual version if it is available
        self.add_software_version(None)

        if refine_json_data_by_sample:
            self._add_general_stats_refine(refine_json_data_by_sample)
            self._add_refine_plot(refine_csv_data_by_sample)

        if cnt_by_cluster_id_by_sample:
            self._add_general_stats_cluster(cnt_by_cluster_id_by_sample)
            self._add_cluster_size_plot(cnt_by_cluster_id_by_sample)

    def _parse_refine(self):
        refine_json_data_by_sample = dict()
        refine_csv_data_by_sample = dict()
        for f in self.find_log_files("isoseq/refine-json", filehandles=True):
            refine_json_data_by_sample[f["s_name"]] = json.load(f["f"])
            self.add_data_source(f, section="refine-json")

        for f in self.find_log_files("isoseq/refine-csv", filehandles=True):
            reader: csv.DictReader = csv.DictReader(f["f"])
            vals_by_metric = defaultdict(list)
            for row in reader:
                for col in ["fivelen", "threelen", "polyAlen", "insertlen"]:
                    vals_by_metric[col].append(int(row[col]))
            refine_csv_data_by_sample[f["s_name"]] = vals_by_metric
            self.add_data_source(f, section="refine-csv")

        refine_json_data_by_sample = self.ignore_samples(refine_json_data_by_sample)
        refine_csv_data_by_sample = self.ignore_samples(refine_csv_data_by_sample)
        log.info(f"Found {len(refine_json_data_by_sample)} refine reports")
        log.info(f"Found {len(refine_csv_data_by_sample)} refine reports")
        if refine_json_data_by_sample:
            self.write_data_file(refine_json_data_by_sample, "multiqc_isoseq_refine_csv")
        if refine_csv_data_by_sample:
            self.write_data_file(refine_csv_data_by_sample, "multiqc_isoseq_refine_json")

        return refine_json_data_by_sample, refine_csv_data_by_sample

    def _parse_cluster(self):
        cnt_by_cluster_id_by_sample = dict()
        for f in self.find_log_files("isoseq/cluster-csv", filehandles=True):
            cnt_by_cluster_id = defaultdict(int)
            reader: csv.DictReader = csv.DictReader(f["f"])
            for row in reader:
                cluster_id = row.get("cluster_id", None)
                if cluster_id is not None:
                    cnt_by_cluster_id[cluster_id] += 1
            if cnt_by_cluster_id:
                cnt_by_cluster_id_by_sample[f["s_name"]] = cnt_by_cluster_id
                self.add_data_source(f)
        if cnt_by_cluster_id_by_sample:
            self.write_data_file(cnt_by_cluster_id_by_sample, "multiqc_isoseq_cluster")
        cnt_by_cluster_id_by_sample = self.ignore_samples(cnt_by_cluster_id_by_sample)
        log.info(f"Found {len(cnt_by_cluster_id_by_sample)} cluster reports")
        if cnt_by_cluster_id_by_sample:
            self.write_data_file(cnt_by_cluster_id_by_sample, "multiqc_isoseq_cluster")
        return cnt_by_cluster_id_by_sample

    def _add_general_stats_cluster(self, size_by_cluster_id_by_sample):
        gstats_data = {}
        for s_name, size_by_cluster_id in size_by_cluster_id_by_sample.items():
            gstats_data[s_name] = {}
            gstats_data[s_name]["n_cluster"] = len(size_by_cluster_id)
            gstats_data[s_name]["mean_cluster_size"] = sum(size_by_cluster_id.values()) / len(size_by_cluster_id)

        headers = {
            "n_cluster": {
                "title": "Clusters",
                "description": "Number of clusters created during the clustering. (1 cluster = 1 transcript)",
                "format": "{:,.d}",
                "scale": "Spectral",
            },
            "mean_cluster_size": {
                "title": "Mean cluster size",
                "description": "Average number of CCS clustered to form transcripts.",
                "format": "{:,.2f}",
                "scale": "RdYlGn",
            },
        }

        self.general_stats_addcols(gstats_data, headers)

    def _add_cluster_size_plot(self, cnt_by_cluster_id_by_sample):
        plot_data = dict()

        for s_name, size_by_cluster_id in cnt_by_cluster_id_by_sample.items():
            plot_data[s_name] = {"=2": 0, "3-10": 0, "11-100": 0, ">100": 0}
            value_counts = defaultdict(int)

            # Calculating value counts similar to df.value_counts("n_CCS") in pandas
            for key, size in size_by_cluster_id.items():
                value_counts[size] += 1

            for n_CCS, size in value_counts.items():
                if n_CCS == 2:
                    plot_data[s_name]["=2"] = size
                elif 2 < n_CCS < 11:
                    plot_data[s_name]["3-10"] += size
                elif 10 < n_CCS < 101:
                    plot_data[s_name]["11-100"] += size
                elif n_CCS > 100:
                    plot_data[s_name][">100"] += size

        cats = {
            "=2": {
                "name": "2 CCS",
            },
            "3-10": {
                "name": "3 to 10 CCS",
            },
            "11-100": {
                "name": "11 to 100 CCS",
            },
            ">100": {
                "name": "More than 100 CCS",
            },
        }

        self.add_section(
            name="Cluster size distribution",
            anchor="isoseq-cluster-size-distribution",
            description="A distribution of cluster size (number of CC clustered to form one Hifi read).",
            helptext="""
            The csv file produced by Iso-Seq cluster shows which CCS have been clustered together to form 
            one Hifi reads. The bargraph represent the distribution of the cluster size using four 
            categories : -2, 3-10, 11-100, >100.
            """,
            plot=bargraph.plot(
                plot_data,
                cats,
                {
                    "id": "isoseq-cluster-size-distribution-barplot",
                    "title": "Iso-Seq cluster: Histogram of cluster size.",
                    "ylab": "Count",  # Y axis label
                },
            ),
        )

        # boxplot_data = {
        #     sname: list(cnt_by_cluster_id.values()) for sname, cnt_by_cluster_id in cnt_by_cluster_id_by_sample.items()
        # }
        # self.add_section(
        #     name="Cluster size distribution",
        #     anchor="isoseq-cluster-size-distribution",
        #     description="A distribution of cluster size (number of CC clustered to form one Hifi read).",
        #     helptext="""
        #     The csv file produced by Iso-Seq cluster shows which CCS have been clustered together to form
        #     one Hifi reads.
        #     """,
        #     plot=box.plot(
        #         boxplot_data,
        #         pconfig={
        #             "id": "isoseq-cluster-size-distribution-barplot",  # HTML ID used for plot
        #             "title": "Iso-Seq cluster: Histogram of cluster size.",
        #         },
        #     ),
        # )

    def _add_general_stats_refine(self, data_by_sample):
        headers = {
            "num_reads_fl": {
                "title": "Full-length",
                "description": "Number of CCS where both primers have been detected",
                "scale": "GnBu",
                "format": "{:,.d}",
            },
            "num_reads_flnc": {
                "title": "Non-chimeric full-length",
                "description": "Number of non-chimeric CCS where both primers have been detected",
                "scale": "RdYlGn",
                "format": "{:,.d}",
            },
            "num_reads_flnc_polya": {
                "title": "Poly(A) free non-chimeric full-length",
                "description": "Number of non-chimeric CCS where both primers have been detected and the poly(A) tail has been removed",
                "scale": "GnBu",
                "format": "{:,.d}",
            },
        }
        self.general_stats_addcols(data_by_sample, headers, namespace="refine")

    def _add_refine_plot(self, values_by_metric_by_sample):
        data_labels = {
            "fivelen": "5' primer length",
            "threelen": "3' primer length",
            "insertlen": "Insert length",
            "polyAlen": "Poly(A) length",
        }
        data_by_sample_by_metric = {k: {} for k in data_labels.keys()}
        for sname, values_by_metric in values_by_metric_by_sample.items():
            for metric in data_labels.keys():
                data_by_sample_by_metric[metric][sname] = values_by_metric[metric]

        data_labels = [
            {
                "label": title,
                "title": "Iso-Seq refine: " + title,
            }
            for metric, title in data_labels.items()
        ]

        self.add_section(
            name="Iso-Seq refine",
            anchor="insert-refine-stats",
            description="Iso-Seq refine statistics",
            helptext="""
            The <code>.report.csv</code> files contain information about 5' prime and 3' primers length, 
            insert length, poly(A) length, and couple of primers detected for each CCS.
            The boxplot presents the distribution based on the known min, max, mean, standard deviation
            statistics for each parameter.
            """,
            plot=box.plot(
                list_of_data_by_sample=list(data_by_sample_by_metric.values()),
                pconfig={
                    "id": "isoseq_refine_boxplot",
                    "title": "Iso-Seq: refine",
                    "data_labels": data_labels,
                },
            ),
        )
