#!/usr/bin/env python

""" MultiQC module to parse output from checkatlas """

from __future__ import print_function

import logging
import os.path
from collections import OrderedDict

from multiqc import config
from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.plots import beeswarm, linegraph, scatter, table

# Initialise the logger
log = logging.getLogger(__name__)

FIG_PATH = "checkatlas_fig"

CELLINDEX_HEADER = "cell_index"
QC_HEADER = ["total_counts", "n_genes_by_counts", "pct_counts_mt"]
QC_RANK_HEADER = ["cellrank_total_counts", "cellrank_n_genes_by_counts", "cellrank_pct_counts_mt"]

LIST_PATTERN = [
    "checkatlas/summary",
    "checkatlas/adata",
    "checkatlas/qc",
    "checkatlas/cluster",
    "checkatlas/annotation",
    "checkatlas/dimred",
    "checkatlas/specificity",
]

DICT_NAMING = {
    "checkatlas/summary": "checkatlas_summ",
    "checkatlas/adata": "_checkatlas_adata",
    "checkatlas/qc": "checkatlas_qc",
    "checkatlas/cluster": "_checkatlas_mcluster",
    "checkatlas/annotation": "checkatlas_mannot",
    "checkatlas/dimred": "_checkatlas_mdimred",
    "checkatlas/specificity": "checkatlas_mspecificity",
}


class MultiqcModule(BaseMultiqcModule):
    """checkatlas module, parses stderr logs."""

    def __init__(self):
        # Initialise the parent object
        super(MultiqcModule, self).__init__(
            name="CheckAtlas",
            anchor="checkatlas",
            href="https://github.com/becavin-lab/checkatlas",
            info="A one-liner tool for quality control of your single-cell atlases.",
            doi="",
        )

        self.data_summary = dict()
        for f in self.find_log_files("checkatlas/summary"):
            s_name = f["s_name"]
            self.data_summary[s_name] = parse_firstline_table_logs(f["f"])
            self.add_data_source(f, s_name)

        self.data_adata = dict()
        for f in self.find_log_files("checkatlas/adata"):
            s_name = f["s_name"]
            self.data_adata[s_name] = parse_firstline_table_logs(f["f"])
            self.add_data_source(f, s_name)

        self.data_qc_counts = dict()
        self.data_qc_genes = dict()
        self.data_qc_mito = dict()
        for f in self.find_log_files("checkatlas/qc"):
            s_name = f["s_name"]
            list_data = parse_qc_logs(f["f"])
            self.data_qc_counts[s_name] = list_data[0]
            self.data_qc_genes[s_name] = list_data[1]
            self.data_qc_mito[s_name] = list_data[2]
            self.add_data_source(f, s_name)

        self.data_metric_cluster = dict()
        for f in self.find_log_files("checkatlas/cluster"):
            s_name = f["s_name"]
            data = parse_metric_logs(f["f"])
            for key, item in data.items():
                self.data_metric_cluster[key] = item
            self.add_data_source(f, s_name)

        self.data_metric_annot = dict()
        for f in self.find_log_files("checkatlas/annotation"):
            s_name = f["s_name"]
            data = parse_metric_logs(f["f"])
            for key, item in data.items():
                self.data_metric_annot[key] = item
            self.add_data_source(f, s_name)

        self.data_metric_dimred = dict()
        for f in self.find_log_files("checkatlas/dimred"):
            s_name = f["s_name"]
            data = parse_metric_logs(f["f"])
            for key, item in data.items():
                self.data_metric_dimred[key] = item
            self.add_data_source(f, s_name)

        # Print results and create sections
        if len(self.data_summary) > 0:
            log.info("Found {} summary tables".format(len(self.data_summary)))
            self.add_summary_section()
        if len(self.data_qc_counts) > 0:
            log.info("Found {} QC counts tables".format(len(self.data_qc_counts)))
            self.add_qc_counts_section()
        if len(self.data_qc_genes) > 0:
            log.info("Found {} QC genes tables".format(len(self.data_qc_genes)))
            self.add_qc_ngenes_section()
        if len(self.data_qc_mito) > 0:
            log.info("Found {} QC mito tables".format(len(self.data_qc_mito)))
            self.add_qc_mito_section()
        if len(self.data_metric_cluster) > 0:
            log.info("Found {} metric cluster tables".format(len(self.data_metric_cluster)))
            self.add_clustermetrics_section()
        if len(self.data_metric_annot) > 0:
            log.info("Found {} metric annot tables".format(len(self.data_metric_annot)))
            self.add_annotationmetrics_section()
        if len(self.data_metric_dimred) > 0:
            log.info("Found {} metric dimred tables".format(len(self.data_metric_dimred)))
            self.add_dimredmetrics_section()
        if len(self.data_adata) > 0:
            log.info("Found {} adata tables".format(len(self.data_adata)))
            self.add_adata_section()

        # Save parsed table
        self.write_data_file(self.data_summary, "multiqc_checkatlas_summary")

    def add_summary_section(self):
        self.add_section(
            name="Atlas overview",
            anchor="checkatlas-summary",
            description="Overview of your single-cell atlases",
            helptext="""
                
                """,
            content=table.plot(self.data_summary),
        )

    def add_adata_section(self):
        config_adata = {"namespace": "adata_table"}
        headers = OrderedDict()
        headers["atlas_obs"] = {"scale": False}
        headers["obsm"] = {"scale": False}
        headers["var"] = {"scale": False}
        headers["varm"] = {"scale": False}
        headers["uns"] = {"scale": False}
        self.add_section(
            name="Atlas object explorer",
            anchor="checkatlas-anndata",
            description="Exploration of your Atlas objects (Scanpy, Cellanger, Seurat)",
            helptext="""
                """,
            content=table.plot(self.data_adata, headers, pconfig=config_adata),
        )

    def add_qc_counts_section(self):
        config_qc = {
            # Building the plot
            "title": "Checkatlas: QC total_counts",
            "ylab": "total_counts",  # X axis label
            "xlab": "log10(Cell Rank)",  # Y axis label
            "id": "qc_counts",  # HTML ID used for plot
            "categories": False,  # Set to True to use x values as categories instead of numbers.
        }

        self.add_section(
            name="QC total_counts",
            anchor="checkatlas-qc_counts",
            description="QC of your atlases log10(Cellrank vs total-counts.",
            helptext="""
            
                """,
            content=linegraph.plot(data=self.data_qc_counts, pconfig=config_qc),
        )

    def add_qc_ngenes_section(self):
        config_qc = {
            # Building the plot
            "title": "Checkatlas: QC n_genes_by_counts",
            "ylab": "n_genes_by_counts",  # X axis label
            "xlab": "log10(Cell Rank)",  # Y axis label
            "logswitch": True,
            "logswitch_active": True,
            "id": "qc_genes",  # HTML ID used for plot
            "categories": False,  # Set to True to use x values as categories instead of numbers.
        }
        self.add_section(
            name="QC n_genes_by_counts",
            anchor="checkatlas-qc_genes",
            description="QC of your atlases log10(Cellrank vs n_genes_by_counts.",
            helptext="""

        """,
            content=linegraph.plot(data=self.data_qc_genes, pconfig=config_qc),
        )

    def add_qc_mito_section(self):
        config_qc = {
            # Building the plot
            "title": "Checkatlas: QC pct_counts_mt",
            "ylab": "pct_counts_mt",  # X axis label
            "xlab": "log10(Cell Rank)",  # Y axis label
            "logswitch": True,
            "logswitch_active": True,
            "id": "qc_mito",  # HTML ID used for plot
            "categories": False,  # Set to True to use x values as categories instead of numbers.
        }
        self.add_section(
            name="QC pct_counts_mt",
            anchor="checkatlas-qc_mito",
            description="QC of your atlases log10(Cellrank vs pct_counts_mt.",
            helptext="",
            content=linegraph.plot(data=self.data_qc_mito, pconfig=config_qc),
        )

    def add_clustermetrics_section(self):
        config_cluster = {"namespace": "metric_cluster_table"}

        self.add_section(
            name="Classification metrics",
            anchor="checkatlas-clustmetrics",
            description="Quality control metrics calculated on your atlases.",
            helptext="""
            
            """,
            content=table.plot(data=self.data_metric_cluster, pconfig=config_cluster),
        )

    def add_annotationmetrics_section(self):
        self.add_section(
            name="Annotation metrics",
            anchor="checkatlas-annotmetrics",
            description="Quality control metrics calculated on your atlases.",
            helptext="""
            
            """,
            content=table.plot(self.data_metric_annot),
        )

    def add_dimredmetrics_section(self):
        self.add_section(
            name="Dimensionality reduction metrics",
            anchor="checkatlas-dimredmetrics",
            description="Quality control metrics calculated on your atlases.",
            helptext="""
                
                """,
            content=table.plot(self.data_metric_dimred),
        )

    def add_specificitymetrics_section(self):
        self.add_section(
            name="Specificity metrics",
            anchor="checkatlas-specmetrics",
            description="Quality control metrics calculated on your atlases.",
            helptext="""
            
            """,
            content="",
        )


def parse_qc_logs(f):
    """
    Parse logs from QC tables in .tsv files
    Order by CellRank
    Calc log10(CelllRank)
    :param f:
    :return:
    """
    lines = f.splitlines()
    headers = lines[0].split("\t")
    index_cellid = headers.index(CELLINDEX_HEADER)
    try:
        index_counts = headers.index(QC_HEADER[0])
        index_rank_counts = headers.index(QC_RANK_HEADER[0])
    except ValueError:
        index_counts = -1
    try:
        index_genes = headers.index(QC_HEADER[1])
        index_rank_genes = headers.index(QC_RANK_HEADER[1])
    except ValueError:
        index_genes = -1
    try:
        index_mito = headers.index(QC_HEADER[2])
        index_rank_mito = headers.index(QC_RANK_HEADER[2])
    except ValueError:
        index_mito = -1

    # get data
    dict_qc_counts = dict()
    dict_qc_genes = dict()
    dict_qc_mito = dict()
    for i in range(1, len(lines)):
        line = lines[i].split("\t")
        cellid = int(line[index_rank_counts])
        if index_counts != -1:
            rank_count = int(line[index_rank_counts])
            count = float(line[index_counts])
            dict_qc_counts[rank_count] = count
        if index_genes != -1:
            rank_genes = float(line[index_rank_genes])
            genes = float(line[index_genes])
            dict_qc_genes[rank_genes] = genes
        if index_mito != -1:
            rank_mito = float(line[index_rank_mito])
            if line[index_mito] != '':
                mito = float(line[index_mito])
            else:
                mito = 0
            dict_qc_mito[rank_mito] = mito

    # reorder by rank
    list_qc_rank_counts = list(dict_qc_counts.keys())
    list_qc_rank_counts.sort()
    list_qc_rank_genes = list(dict_qc_genes.keys())
    list_qc_rank_genes.sort()
    list_qc_rank_mito = list(dict_qc_mito.keys())
    list_qc_rank_mito.sort()
    data_qc_counts = dict()
    data_qc_genes = dict()
    data_qc_mito = dict()
    for rank in list_qc_rank_counts:
        data_qc_counts[rank] = dict_qc_counts[rank]
    for rank in list_qc_rank_genes:
        data_qc_genes[rank] = dict_qc_genes[rank]
    for rank in list_qc_rank_mito:
        data_qc_mito[rank] = dict_qc_mito[rank]
    list_data = [
        data_qc_counts,
        data_qc_genes,
        data_qc_mito,
    ]
    return list_data


def parse_firstline_table_logs(f):
    """
    Parse only header and first line of logs which are .tsv files
    Ex: _checkatlas.tsv, adata_checkatlas.tsv
    :param f:
    :return:
    """
    data = {}
    lines = f.splitlines()
    headers = lines[0].split("\t")
    for i in range(1, len(lines)):
        line = lines[i].split("\t")
        for j in range(0, len(line)):
            data[headers[j]] = line[j]
    return data


def parse_metric_logs(f):
    """
    Parse all lines of logs which are .tsv files
    Ex: _checkatlas.tsv, adata_checkatlas.tsv
    :param f:
    :return:
    """
    data = {}
    lines = f.splitlines()
    headers = lines[0].split("\t")
    for i in range(1, len(lines)):
        line = lines[i].split("\t")
        line_dict = {}
        for j in range(1, len(line)):
            line_dict[headers[j]] = line[j]
        data[line[0]] = line_dict
    return data
