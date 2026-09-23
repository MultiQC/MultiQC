"""`simplex` / `duplex` / `codec --stats`: key/value/description rows of consensus-calling statistics."""

import logging
from typing import Any, Dict, Set

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

from .schemas import ConsensusStatMetric
from .util import drop_none, load_rows, pct, sample_name

log = logging.getLogger(__name__)

REJECTED_PREFIX = "raw_reads_rejected_for_"
# Keys every fgumi consensus caller writes, read directly for the General Statistics columns.
REQUIRED_KEYS = ("consensus_reads_emitted", "frac_raw_reads_used")


_GENERAL_STATS_HEADERS: Dict[str, Dict[str, Any]] = {
    "consensus_reads_emitted": {
        "title": "Consensus reads",
        "description": "Consensus reads emitted",
        "format": "{:,.0f}",
        "scale": "Greens",
    },
    "frac_raw_reads_used": {
        "title": "% reads used",
        "description": "Raw reads used in a consensus read",
        "suffix": "%",
        "max": 100,
        "min": 0,
        "hidden": True,
    },
}


def parse_reports(module: BaseMultiqcModule) -> Set[str]:
    stats: Dict[str, Dict[str, float]] = {}
    for f in module.find_log_files("fgumi/consensus_stats"):
        rows = load_rows(f, ConsensusStatMetric)
        if rows is None:
            continue
        values: Dict[str, float] = {}
        for row in rows:
            try:
                values[row.key] = float(row.value)
            except ValueError:
                continue
        if "raw_reads_considered" not in values:
            log.debug(f"{f['fn']}: key/value file is not an fgumi consensus --stats file, skipping")
            continue
        missing = [key for key in REQUIRED_KEYS if key not in values]
        if missing:
            log.warning(f"Skipping {f['fn']}: fgumi consensus --stats file is missing {missing}")
            continue
        s_name = sample_name(module, f)
        module.add_data_source(f, s_name)
        module.add_software_version(None, s_name)
        stats[s_name] = values
    stats = module.ignore_samples(stats)
    if not stats:
        return set()

    bars = {
        s: {k: v for k, v in values.items() if (k == "raw_reads_used" or k.startswith(REJECTED_PREFIX)) and v > 0}
        for s, values in stats.items()
    }
    cats = {"raw_reads_used": {"name": "Used in a consensus"}}
    for key in sorted({k for values in bars.values() for k in values if k.startswith(REJECTED_PREFIX)}):
        cats[key] = {"name": "Rejected: " + key[len(REJECTED_PREFIX) :].replace("_", " ")}
    module.add_section(
        name="Consensus calling",
        anchor="fgumi-consensus-stats",
        description="Raw reads used in consensus reads, and why the rest were rejected.",
        helptext="Written by `fgumi simplex`, `duplex` or `codec` with `--stats`.",
        plot=bargraph.plot(
            bars, cats, {"id": "fgumi_consensus_rejections", "title": "fgumi: Consensus calling", "ylab": "Raw reads"}
        ),
    )
    module.general_stats_addcols(
        {
            s: drop_none(
                {
                    "consensus_reads_emitted": v["consensus_reads_emitted"],
                    "frac_raw_reads_used": pct(v["frac_raw_reads_used"]),
                }
            )
            for s, v in stats.items()
        },
        _GENERAL_STATS_HEADERS,
    )
    module.write_data_file(stats, "multiqc_fgumi_consensus_stats")
    return set(stats)
