"""`review`: the `<output>.txt` detail file of consensus reads supporting each reviewed variant, summarized."""

from typing import Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import table

from .schemas import MetricFormatError, ReviewDetailMetric
from .util import register, skip_unreadable, stream_dicts


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    data: Dict[str, Dict[str, int]] = {}
    for f in module.find_log_files("fgumi/review", filehandles=True):
        if f["f"] is None:
            continue
        sites, reads, observations = set(), set(), 0
        try:
            for row in stream_dicts(f["f"], f["fn"], ReviewDetailMetric):
                sites.add((row["chrom"], row["pos"]))
                reads.add(row["consensus_read"])
                observations += 1
        except (MetricFormatError, ValueError) as error:  # UnicodeDecodeError is a ValueError
            skip_unreadable(f, error)
            continue
        data[register(module, f)] = {"sites": len(sites), "consensus_reads": len(reads), "observations": observations}
    data = module.ignore_samples(data)
    if not data:
        return set()
    module.add_section(
        name="Consensus variant review",
        anchor="fgumi-review",
        description="Consensus reads with a non-reference base at a variant reviewed by `fgumi review` (or fgbio "
        "ReviewConsensusVariants).",
        helptext="Summarized from the `<output>.txt` detail file, which has a row only for each consensus read with a "
        "non-reference, non-deleted base at a reviewed site. Sites where no consensus read differs from the reference "
        "are not in the file, so the number of sites reviewed cannot be recovered from it. Consensus reads count R1 "
        "and R2 separately.",
        plot=table.plot(
            data,
            {
                "sites": {
                    "title": "Sites with non-ref reads",
                    "description": "Reviewed sites with at least one non-reference consensus read",
                    "format": "{:,.0f}",
                    "scale": "Blues",
                },
                "consensus_reads": {
                    "title": "Non-ref consensus reads",
                    "description": "Distinct consensus read ends (R1 and R2 counted separately) with a non-reference base",
                    "format": "{:,.0f}",
                    "scale": "Greens",
                },
                "observations": {
                    "title": "Non-ref observations",
                    "description": "Rows of the detail file: one per consensus read and site",
                    "format": "{:,.0f}",
                    "scale": "Purples",
                },
            },
            {"id": "fgumi_review", "title": "fgumi: Consensus variant review", "namespace": "fgumi review"},
        ),
    )
    module.write_data_file(data, "multiqc_fgumi_review")
    return set(data)
