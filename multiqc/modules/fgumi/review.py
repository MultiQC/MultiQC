"""`review`: the `<output>.txt` detail file of consensus reads supporting each reviewed variant, summarized."""

import logging
from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

from .schemas import MetricFormatError, ReviewDetailMetric
from .util import sample_name, stream_dicts

log = logging.getLogger(__name__)


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, int]] = {}
    for f in module.find_log_files("fgumi/review", filehandles=True):
        sites, reads, observations = set(), set(), 0
        try:
            for row in stream_dicts(f, ReviewDetailMetric):
                sites.add((row["chrom"], row["pos"]))
                reads.add(row["consensus_read"])
                observations += 1
        except MetricFormatError as error:
            log.warning(f"Skipping {error}")
            continue
        s_name = sample_name(module, f)
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        data[s_name] = {"sites": len(sites), "consensus_reads": len(reads), "observations": observations}
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Consensus variant review",
        anchor="fgumi-review",
        description="Variant sites reviewed by `fgumi review`, and the consensus reads observed at them.",
        helptext="Summarized from the `<output>.txt` detail file of `fgumi review`, which has one row per variant site and "
        "consensus read.",
        plot=table.plot(
            data,
            {
                "sites": {"title": "Variant sites", "format": "{:,.0f}", "scale": "Blues"},
                "consensus_reads": {"title": "Consensus reads", "format": "{:,.0f}", "scale": "Greens"},
                "observations": {"title": "Observations", "format": "{:,.0f}", "scale": "Purples"},
            },
            {"id": "fgumi_review", "title": "fgumi: Consensus variant review", "namespace": "fgumi review"},
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_review")
    return set(data)
