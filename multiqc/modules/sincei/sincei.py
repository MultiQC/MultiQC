import logging

from multiqc.base_module import BaseMultiqcModule, ModuleNoSamplesFound

# sincei modules
from .scCountQC import scCountQCMixin
from .scFilterStats import scFilterStatsMixin

# Initialise the logger
log = logging.getLogger(__name__)


class MultiqcModule(
    BaseMultiqcModule,
    scFilterStatsMixin,
    scCountQCMixin,
):
    """
    sincei (short for Single Cell Informatics) is a command-line toolkit for
    exploration of single-cell epigenomics data. It accommodates data from a
    wide range of single-cell protocols, such as droplet-based (10x Genomics)
    and plate-based protocols, gene expression (scRNA-seq) and epigenomics
    (scATAC, scCUTnTAG, scBS-seq). sincei can be used for quality control of
    these datasets directly from BAM files (read-level QC), as well as after
    signal aggregation (count-level). Additional functionalities include
    filtering, dimensionality reduction, and plotting of single-cell data.

    The MultiQC module for sincei parses the following text outputs:

    - `scFilterStats` (default output file)
    - `scCountQC --outMetrics` (currently the cell-level metrics are supported)

    sincei reports one row per cell, identified by `Cell_ID`. MultiQC parses
    the sample name from the part of `Cell_ID` before `::` and aggregates
    metrics across cells by taking the median per sample. Each row in the
    report tables therefore represents one sample, summarising all of its
    cells.
    """

    def __init__(self):
        # Initialise the parent object
        super().__init__(
            name="sincei",
            anchor="sincei",
            target="sincei",
            href="https://sincei.readthedocs.io",
            info="Toolkit for processing and analyzing single-cell (epi)genomics data.",
            doi="10.5281/zenodo.7853375",
        )

        samples = dict()
        cells = dict()

        # scFilterStats
        samples["scFilterStats"], cells["scFilterStats"] = self.parse_scFilterStats()
        if samples["scFilterStats"] > 0:
            log.debug(f"Found {samples['scFilterStats']} sincei scFilterStats samples ({cells['scFilterStats']} cells)")

        # scCountQC
        samples["scCountQC"], cells["scCountQC"] = self.parse_scCountQC()
        if samples["scCountQC"] > 0:
            log.debug(f"Found {samples['scCountQC']} sincei scCountQC samples ({cells['scCountQC']} cells)")

        total_samples = sum(samples.values())
        total_cells = sum(cells.values())
        if total_samples > 0:
            log.info(f"Found {total_samples} samples ({total_cells} cells)")
        else:
            raise ModuleNoSamplesFound
