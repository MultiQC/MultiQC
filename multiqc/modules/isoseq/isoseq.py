# !/usr/bin/env python

"""
MultiQC module to parse output from Iso-Seq tools
"""

import logging

from multiqc.modules.base_module import BaseMultiqcModule
from multiqc.modules.isoseq.cluster import ClusterMixin
from multiqc.modules.isoseq.refine import RefineMixin

log = logging.getLogger(__name__)


class MultiqcModule(BaseMultiqcModule, RefineMixin, ClusterMixin):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name="Iso-Seq",
            anchor="isoseq",
            href="https://github.com/PacificBiosciences/IsoSeq",
            info="contains the newest tools to identify transcripts in PacBio single-molecule sequencing data (HiFi reads).",
            # doi=,  # Not published
        )

        refine_data = self.parse_refine_log()
        refine_data = self.ignore_samples(refine_data)
        if refine_data:
            log.info(f"Found {len(refine_data)} refine reports")

        cluster_data = self.parse_cluster_log()
        cluster_data = self.ignore_samples(cluster_data)
        if cluster_data:
            log.info(f"Found {len(cluster_data)} cluster reports")

        # If we found no data
        if not refine_data and not cluster_data:
            raise UserWarning

        if refine_data:
            self.add_general_stats_refine(refine_data)
            self.add_section_fivelen(refine_data)
            self.add_section_insertlen(refine_data)
            self.add_section_polyalen(refine_data)
            self.add_section_threelen(refine_data)

        if cluster_data:
            self.add_general_stats_cluster(cluster_data)
            self.add_section_cluster_size_bargraph(cluster_data)
