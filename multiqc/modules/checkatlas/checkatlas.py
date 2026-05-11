import logging
from typing import Union

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound
from multiqc.plots import linegraph, table

log = logging.getLogger(__name__)


# Maps each per-cell QC value column to its companion "cell rank" column in checkatlas qc/*.tsv.
QC_METRICS: dict[str, str] = {
    "total_counts": "cellrank_total_counts",
    "n_genes_by_counts": "cellrank_n_genes_by_counts",
    "pct_counts_mt": "cellrank_pct_counts_mt",
}

# Cap the number of points plotted per atlas in the QC line graphs.
# Real atlases can have hundreds of thousands of cells; rendering one trace per cell
# makes the report sluggish without adding visual information.
QC_PLOT_MAX_POINTS = 1000

Coerced = Union[int, float, str]


class MultiqcModule(BaseMultiqcModule):
    """
    CheckAtlas is a one-liner tool to check the quality of your single-cell atlases. For every atlas, it produces
    quality control tables and figures which can then be processed by MultiQC. CheckAtlas is able to load Scanpy,
    Seurat, and CellRanger files.

    MultiQC parses the following tables produced by CheckAtlas:

    - `summary/sample_name.tsv` - Summary tables with general information on atlases
    - `adata/sample_name.tsv` - Table with all scanpy adata features
    - `qc/sample_name.tsv` - Quality control tables for every atlas
    - `cluster/sample_name.tsv` - Table with cluster metrics calculated for every atlas
    - `annot/sample_name.tsv` - Table with annotation metrics calculated for every atlas
    - `dimred/sample_name.tsv` - Table with dimensionality reduction metrics calculated for every atlas
    """

    def __init__(self):
        super().__init__(
            name="CheckAtlas",
            anchor="checkatlas",
            href="https://github.com/becavin-lab/checkatlas",
            info="A one-liner tool for quality control of your single-cell atlases.",
            # Can't find a DOI // doi=
        )

        self.data_summary: dict[str, dict[str, Coerced]] = {}
        for f in self.find_log_files("checkatlas/summary"):
            s_name = f["s_name"]
            self.data_summary[s_name] = parse_first_row(f["f"])
            self.add_data_source(f, s_name)

        self.data_adata: dict[str, dict[str, Coerced]] = {}
        for f in self.find_log_files("checkatlas/adata"):
            s_name = f["s_name"]
            self.data_adata[s_name] = parse_first_row(f["f"])
            self.add_data_source(f, s_name)

        self.data_qc_counts: dict[str, dict[int, float]] = {}
        self.data_qc_genes: dict[str, dict[int, float]] = {}
        self.data_qc_mito: dict[str, dict[int, float]] = {}
        for f in self.find_log_files("checkatlas/qc"):
            s_name = f["s_name"]
            parsed = parse_qc_logs(f["f"])
            if parsed["total_counts"]:
                self.data_qc_counts[s_name] = _downsample(parsed["total_counts"])
            if parsed["n_genes_by_counts"]:
                self.data_qc_genes[s_name] = _downsample(parsed["n_genes_by_counts"])
            if parsed["pct_counts_mt"]:
                self.data_qc_mito[s_name] = _downsample(parsed["pct_counts_mt"])
            self.add_data_source(f, s_name)

        self.data_metric_cluster = self._load_metric_table("checkatlas/cluster")
        self.data_metric_annot = self._load_metric_table("checkatlas/annotation")
        self.data_metric_dimred = self._load_metric_table("checkatlas/dimred")

        # Filter ignored samples from every dataset before rendering anything.
        self.data_summary = self.ignore_samples(self.data_summary)
        self.data_adata = self.ignore_samples(self.data_adata)
        self.data_qc_counts = self.ignore_samples(self.data_qc_counts)
        self.data_qc_genes = self.ignore_samples(self.data_qc_genes)
        self.data_qc_mito = self.ignore_samples(self.data_qc_mito)
        self.data_metric_cluster = self.ignore_samples(self.data_metric_cluster)
        self.data_metric_annot = self.ignore_samples(self.data_metric_annot)
        self.data_metric_dimred = self.ignore_samples(self.data_metric_dimred)

        all_data: list[dict] = [
            self.data_summary,
            self.data_adata,
            self.data_qc_counts,
            self.data_qc_genes,
            self.data_qc_mito,
            self.data_metric_cluster,
            self.data_metric_annot,
            self.data_metric_dimred,
        ]
        if not any(all_data):
            raise ModuleNoSamplesFound

        self.add_software_version(None)

        if self.data_summary:
            log.info(f"Found {len(self.data_summary)} summary tables")
            self.add_summary_section()
        if self.data_qc_counts:
            log.info(f"Found {len(self.data_qc_counts)} QC counts tables")
            self.add_qc_counts_section()
        if self.data_qc_genes:
            log.info(f"Found {len(self.data_qc_genes)} QC genes tables")
            self.add_qc_ngenes_section()
        if self.data_qc_mito:
            log.info(f"Found {len(self.data_qc_mito)} QC mito tables")
            self.add_qc_mito_section()
        if self.data_metric_cluster:
            log.info(f"Found {len(self.data_metric_cluster)} cluster metric tables")
            self.add_clustermetrics_section()
        if self.data_metric_annot:
            log.info(f"Found {len(self.data_metric_annot)} annotation metric tables")
            self.add_annotationmetrics_section()
        if self.data_metric_dimred:
            log.info(f"Found {len(self.data_metric_dimred)} dimred metric tables")
            self.add_dimredmetrics_section()
        if self.data_adata:
            log.info(f"Found {len(self.data_adata)} adata tables")
            self.add_adata_section()

        # Write one data file per non-empty section.
        for fname, data in [
            ("multiqc_checkatlas_summary", self.data_summary),
            ("multiqc_checkatlas_adata", self.data_adata),
            ("multiqc_checkatlas_qc_counts", self.data_qc_counts),
            ("multiqc_checkatlas_qc_genes", self.data_qc_genes),
            ("multiqc_checkatlas_qc_mito", self.data_qc_mito),
            ("multiqc_checkatlas_cluster", self.data_metric_cluster),
            ("multiqc_checkatlas_annotation", self.data_metric_annot),
            ("multiqc_checkatlas_dimred", self.data_metric_dimred),
        ]:
            if data:
                self.write_data_file(data, fname)

    def _load_metric_table(self, search_key: str) -> dict[str, dict[str, Coerced]]:
        """Load all metric .tsv files for a given search key, warning on key collisions."""
        data: dict[str, dict[str, Coerced]] = {}
        for f in self.find_log_files(search_key):
            parsed = parse_metric_logs(f["f"])
            duplicate_keys = [key for key in parsed if key in data]
            data.update(parsed)
            if duplicate_keys:
                first = duplicate_keys[0]
                extra = f" (and {len(duplicate_keys) - 1} more)" if len(duplicate_keys) > 1 else ""
                log.warning(
                    f"{f['fn']}: {len(duplicate_keys)} duplicate {search_key} key(s) overwrote earlier values; "
                    f"first was '{first}'{extra}"
                )
            self.add_data_source(f, f["s_name"])
        return data

    def add_summary_section(self):
        pconfig = {
            "namespace": "summary_table",
            "id": "checkatlas_summary",
            "title": "CheckAtlas: Atlas Summary",
        }
        headers = {
            "AtlasFileType": {
                "title": "Atlas Type",
                "description": "Source file format of the atlas (e.g. AnnData, Seurat, CellRanger)",
            },
            "NbCells": {
                "title": "Cells",
                "description": "Number of cells in the atlas",
                "format": "{:,.0f}",
                "scale": "Blues",
            },
            "NbGenes": {
                "title": "Genes",
                "description": "Number of genes (features) in the atlas",
                "format": "{:,.0f}",
                "scale": "Greens",
            },
            "AnnData.raw": {
                "title": "Has Raw Layer",
                "description": "Whether the atlas contains a raw counts layer (adata.raw)",
            },
            "AnnData.X": {
                "title": "Has X Layer",
                "description": "Whether the atlas contains the main expression matrix (adata.X)",
            },
            "File_extension": {
                "title": "Atlas File",
                "description": (
                    "Source file or extension for the atlas "
                    "(file name for AnnData, extension for Seurat, directory name for CellRanger)"
                ),
            },
            "File_path": {
                "title": "File Path",
                "description": "Path to the source atlas file on disk",
                "hidden": True,
            },
        }
        self.add_section(
            name="Atlas Overview",
            anchor="checkatlas_summary_section",
            description="Top-level statistics for each input atlas.",
            helptext="""
                One row per atlas, showing the source format (AnnData, Seurat or CellRanger),
                the number of cells and genes, and which expression layers (`adata.X`, `adata.raw`)
                are present. Useful for comparing scale and shape of atlases at a glance.
                """,
            plot=table.plot(self.data_summary, headers, pconfig=pconfig),
        )

    def add_adata_section(self):
        pconfig = {
            "namespace": "adata_table",
            "id": "checkatlas_adata",
            "title": "CheckAtlas: Atlas Object Contents",
        }
        headers = {
            "atlas_obs": {
                "title": "Observations (obs)",
                "description": "Per-cell annotation keys stored in adata.obs (e.g. cluster labels, QC metrics)",
            },
            "obsm": {
                "title": "Multi-dim Observations (obsm)",
                "description": "Per-cell multi-dimensional arrays stored in adata.obsm (e.g. PCA, UMAP embeddings)",
            },
            "var": {
                "title": "Variables (var)",
                "description": "Per-gene annotation keys stored in adata.var (e.g. gene IDs, feature types)",
            },
            "varm": {
                "title": "Multi-dim Variables (varm)",
                "description": "Per-gene multi-dimensional arrays stored in adata.varm (e.g. PCA loadings)",
            },
            "uns": {
                "title": "Unstructured (uns)",
                "description": (
                    "Unstructured annotation keys stored in adata.uns "
                    "(e.g. colour maps, clustering parameters, rank-genes results)"
                ),
            },
        }
        self.add_section(
            name="Atlas Object Contents",
            anchor="checkatlas_adata_section",
            description="Annotation keys and slots present inside each atlas object.",
            helptext="""
                For each atlas, the cells of this table list the keys stored in the corresponding
                AnnData slot (Seurat and CellRanger atlases are mapped onto the same structure).
                Use this to confirm that expected metadata (cluster assignments, embeddings,
                QC metrics) is present in each atlas.
                """,
            plot=table.plot(self.data_adata, headers, pconfig=pconfig),
        )

    def add_qc_counts_section(self):
        pconfig = {
            "title": "CheckAtlas: Total Counts per Cell",
            "ylab": "Total Counts per Cell",
            "xlab": "Cell Rank",
            "id": "checkatlas_qc_counts",
            "categories": False,
        }
        self.add_section(
            name="QC: Counts per Cell",
            anchor="checkatlas_qc_counts_section",
            description="Total counts per cell, ordered from highest to lowest.",
            helptext="""
                For each atlas, cells are ranked by their `total_counts` (the sum of UMI or read
                counts across all genes) and plotted from highest to lowest. Comparing the shape
                and height of the curves between atlases highlights differences in sequencing
                depth and capture efficiency.
                """,
            plot=linegraph.plot(self.data_qc_counts, pconfig=pconfig),
        )

    def add_qc_ngenes_section(self):
        pconfig = {
            "title": "CheckAtlas: Genes Detected per Cell",
            "ylab": "Genes Detected per Cell",
            "xlab": "Cell Rank",
            "logswitch": True,
            "logswitch_active": True,
            "id": "checkatlas_qc_genes",
            "categories": False,
        }
        self.add_section(
            name="QC: Genes Detected per Cell",
            anchor="checkatlas_qc_genes_section",
            description="Number of genes detected per cell, ordered from highest to lowest.",
            helptext="""
                For each atlas, cells are ranked by their `n_genes_by_counts` (the number of genes
                with at least one count) and plotted from highest to lowest. Low gene complexity
                may indicate poor capture or low-quality cells; large differences between atlases
                can point to batch effects or platform differences.
                """,
            plot=linegraph.plot(self.data_qc_genes, pconfig=pconfig),
        )

    def add_qc_mito_section(self):
        pconfig = {
            "title": "CheckAtlas: Mitochondrial Counts per Cell",
            "ylab": "Mitochondrial Counts (%)",
            "xlab": "Cell Rank",
            "logswitch": True,
            "logswitch_active": True,
            "id": "checkatlas_qc_mito",
            "categories": False,
        }
        self.add_section(
            name="QC: Mitochondrial Counts",
            anchor="checkatlas_qc_mito_section",
            description="Percentage of counts from mitochondrial genes per cell, ordered from highest to lowest.",
            helptext="""
                For each atlas, cells are ranked by their `pct_counts_mt` (percentage of total
                counts that come from mitochondrial genes) and plotted from highest to lowest.
                High mitochondrial content is a common marker of stressed or dying cells and is
                a standard QC filter for single-cell datasets. This column is only present in
                atlases where mitochondrial genes have been annotated upstream; cells with an
                empty value are skipped.
                """,
            plot=linegraph.plot(self.data_qc_mito, pconfig=pconfig),
        )

    def add_clustermetrics_section(self):
        pconfig = {
            "namespace": "metric_cluster_table",
            "id": "checkatlas_cluster",
            "title": "CheckAtlas: Clustering Metrics",
        }
        headers = {
            "obs": {
                "title": "Clustering Key",
                "description": "Name of the obs column in the atlas holding the cluster assignments",
            },
            "davies_bouldin": {
                "title": "Davies-Bouldin",
                "description": (
                    "Davies-Bouldin index: average similarity of each cluster with its most similar cluster, "
                    "where similarity is the ratio of within-cluster to between-cluster distances. "
                    "Lower is better; 0 is the best possible score."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn-rev",
                "min": 0,
            },
            "silhouette": {
                "title": "Silhouette",
                "description": (
                    "Silhouette score: measures how similar each cell is to its own cluster compared to "
                    "other clusters. Ranges from -1 (poor) to +1 (well clustered); higher is better."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn",
                "min": -1,
                "max": 1,
            },
            "calinski_harabasz": {
                "title": "Calinski-Harabasz",
                "description": (
                    "Calinski-Harabasz index (variance ratio criterion): ratio of between-cluster to "
                    "within-cluster dispersion. Higher values indicate better-defined clusters."
                ),
                "format": "{:,.2f}",
                "scale": "RdYlGn",
                "min": 0,
            },
        }
        self.add_section(
            name="Clustering Metrics",
            anchor="checkatlas_cluster_section",
            description="Clustering quality metrics calculated for each atlas.",
            helptext="""
                Internal clustering evaluation scores (Davies-Bouldin, Silhouette, Calinski-Harabasz)
                computed by CheckAtlas for the clusterings stored in each atlas. Higher-quality
                clusterings generally show better separation between clusters and tighter compactness
                within them.
                """,
            plot=table.plot(self.data_metric_cluster, headers, pconfig=pconfig),
        )

    def add_annotationmetrics_section(self):
        pconfig = {
            "namespace": "metric_annot_table",
            "id": "checkatlas_annot",
            "title": "CheckAtlas: Annotation Metrics",
        }
        zero_to_one = {"format": "{:,.3f}", "scale": "RdYlGn", "min": 0, "max": 1}
        headers = {
            "Reference": {
                "title": "Reference",
                "description": "Reference clustering or label set the annotation is being compared against",
            },
            "obs": {
                "title": "Annotation Key",
                "description": "Name of the obs column in the atlas holding the cell annotations",
            },
            "rand_index": {
                "title": "Rand Index",
                "description": (
                    "Rand index: fraction of cell pairs whose annotation agrees with the reference. "
                    "Ranges from 0 (no agreement) to 1 (identical assignments)."
                ),
                **zero_to_one,
            },
            "adj_rand_index": {
                "title": "Adjusted Rand Index",
                "description": (
                    "Rand index corrected for chance. Ranges from ~0 (random labelling) to 1 "
                    "(perfect agreement); can be slightly negative for worse-than-random labelling."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn",
                "max": 1,
            },
            "mutual_info": {
                "title": "Mutual Information",
                "description": (
                    "Mutual information between the annotation and the reference labels. "
                    "Higher values indicate more shared information."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn",
                "min": 0,
            },
            "adj_mutual_info": {
                "title": "Adjusted Mutual Information",
                "description": "Mutual information corrected for chance. 0 ≈ random labelling, 1 = perfect agreement.",
                **zero_to_one,
            },
            "normalized_mutual_info": {
                "title": "Normalised Mutual Information",
                "description": (
                    "Mutual information normalised to [0, 1]. 0 = no shared information, 1 = identical labellings."
                ),
                **zero_to_one,
            },
            "fowlkes_mallows": {
                "title": "Fowlkes-Mallows",
                "description": (
                    "Geometric mean of pairwise precision and recall between the annotation and the "
                    "reference. Ranges from 0 to 1; higher is better."
                ),
                **zero_to_one,
            },
            "vmeasure": {
                "title": "V-measure",
                "description": (
                    "Harmonic mean of homogeneity and completeness of the annotation with respect to the "
                    "reference. Ranges from 0 to 1; higher is better."
                ),
                **zero_to_one,
            },
        }
        self.add_section(
            name="Annotation Metrics",
            anchor="checkatlas_annot_section",
            description="Cell annotation quality metrics calculated for each atlas.",
            helptext="""
                Metrics evaluating cell-type annotations stored in each atlas, comparing them against
                a reference clustering or label set. CheckAtlas can produce a range of comparison
                metrics (Rand Index, Adjusted Rand Index, Mutual Information variants, Fowlkes-Mallows,
                V-measure, etc.); only the metrics that were computed will appear as columns here.
                """,
            plot=table.plot(self.data_metric_annot, headers, pconfig=pconfig),
        )

    def add_dimredmetrics_section(self):
        pconfig = {
            "namespace": "metric_dimred_table",
            "id": "checkatlas_dimred",
            "title": "CheckAtlas: Dimensionality Reduction Metrics",
        }
        headers = {
            "obsm": {
                "title": "Embedding Key",
                "description": "Name of the obsm slot holding the dimensionality reduction (e.g. X_pca, X_umap)",
            },
            "kruskal_stress": {
                "title": "Kruskal Stress",
                "description": (
                    "Kruskal stress: measures how well pairwise distances in the low-dimensional "
                    "embedding match those in the original high-dimensional space. Lower is better."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn-rev",
                "min": 0,
            },
            "spearman_rho": {
                "title": "Spearman ρ",
                "description": (
                    "Spearman rank correlation between pairwise distances in the high- and "
                    "low-dimensional spaces. Closer to 1 means distances are well preserved."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn",
                "min": -1,
                "max": 1,
            },
            "entourage": {
                "title": "Entourage",
                "description": (
                    "Entourage score: proportion of each cell's nearest neighbours in the high-dimensional "
                    "space that remain its neighbours in the low-dimensional embedding. Higher is better."
                ),
                "format": "{:,.3f}",
                "scale": "RdYlGn",
                "min": 0,
                "max": 1,
            },
        }
        self.add_section(
            name="Dimensionality Reduction Metrics",
            anchor="checkatlas_dimred_section",
            description="Dimensionality reduction quality metrics calculated for each atlas.",
            helptext="""
                Metrics evaluating how well dimensionality reductions (PCA, UMAP, t-SNE, etc.) preserve
                the structure of the original high-dimensional data. CheckAtlas can report Kruskal
                stress, Spearman ρ on pairwise distances, and an Entourage neighbour-preservation
                score; only the metrics that were computed will appear as columns here.
                """,
            plot=table.plot(self.data_metric_dimred, headers, pconfig=pconfig),
        )


def _coerce(value: str) -> Coerced:
    """Convert a TSV field to int / float when possible, otherwise return the original string."""
    if value == "":
        return value
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def parse_qc_logs(f: str) -> dict[str, dict[int, float]]:
    """Parse a QC .tsv table into one rank→value dict per QC metric."""
    lines = f.splitlines()
    headers = lines[0].split("\t")

    indices: dict[str, tuple[int, int]] = {}
    for value_col, rank_col in QC_METRICS.items():
        try:
            indices[value_col] = (headers.index(value_col), headers.index(rank_col))
        except ValueError:
            continue

    series: dict[str, dict[int, float]] = {value_col: {} for value_col in QC_METRICS}
    # start=2 so line_no matches the file line number (line 1 is the header).
    for line_no, raw_line in enumerate(lines[1:], start=2):
        line = raw_line.split("\t")
        for value_col, (value_idx, rank_idx) in indices.items():
            if max(value_idx, rank_idx) >= len(line):
                continue
            value_str = line[value_idx]
            if value_str == "":
                # Skip empty values consistently across metrics.
                log.warning(f"Empty {value_col} value on line {line_no}; skipping")
                continue
            try:
                series[value_col][int(line[rank_idx])] = float(value_str)
            except ValueError:
                log.warning(f"Could not parse {value_col} on line {line_no}: {raw_line!r}")

    return {value_col: dict(sorted(rows.items())) for value_col, rows in series.items()}


def parse_first_row(f: str) -> dict[str, Coerced]:
    """Parse the header and first data row of a .tsv table into a flat dict."""
    lines = f.splitlines()
    if len(lines) < 2:
        return {}
    headers = lines[0].split("\t")
    line = lines[1].split("\t")
    return {headers[j]: _coerce(line[j]) for j in range(min(len(line), len(headers)))}


def parse_metric_logs(f: str) -> dict[str, dict[str, Coerced]]:
    """Parse a .tsv metric table keyed by the first column of each row."""
    data: dict[str, dict[str, Coerced]] = {}
    lines = f.splitlines()
    headers = lines[0].split("\t")
    for raw_line in lines[1:]:
        line = raw_line.split("\t")
        if not line or not line[0]:
            continue
        row = {headers[j]: _coerce(line[j]) for j in range(1, min(len(line), len(headers)))}
        data[line[0]] = row
    return data


def _downsample(series: dict[int, float], max_points: int = QC_PLOT_MAX_POINTS) -> dict[int, float]:
    """Subsample a rank-keyed series down to at most max_points evenly-spaced points."""
    if len(series) <= max_points:
        return series
    items = sorted(series.items())
    step = len(items) / max_points
    sampled = [items[int(i * step)] for i in range(max_points)]
    if sampled[-1] != items[-1]:
        sampled.append(items[-1])
    return dict(sampled)
