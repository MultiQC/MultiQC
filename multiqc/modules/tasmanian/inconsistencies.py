"""MultiQC submodule to parse read_pair_inconsistencies.tsv from tasmanian-diagnostics"""

import logging
import re
from typing import Dict

from multiqc import BaseMultiqcModule

from .util import add_category_bargraph

log = logging.getLogger(__name__)

DISCORDANCE_RE = re.compile(r"^R1:([ACGTN])_R2:([ACGTN])$")


def parse_tasmanian_inconsistencies(module: BaseMultiqcModule) -> int:
    """Find tasmanian-diagnostics read pair inconsistency tables and parse their data"""

    inconsistency_data: Dict[str, Dict] = {}

    for f in module.find_log_files("tasmanian/inconsistencies"):
        try:
            data = parse_inconsistencies_table(f["f"])
        except ValueError as e:
            raise ValueError(f"Could not parse tasmanian-diagnostics output {f['fn']}: {e}") from e

        s_name = module.clean_s_name(f["s_name"], f)
        if s_name in inconsistency_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        module.add_data_source(f, s_name=s_name, section="inconsistencies")
        inconsistency_data[s_name] = data

        # tasmanian-diagnostics does not record its version in the output
        module.add_software_version(None, s_name)

    inconsistency_data = module.ignore_samples(inconsistency_data)
    if len(inconsistency_data) == 0:
        return 0

    headers: Dict = {
        "inconsistent_bases": {
            "title": "Overlap Inconsistencies",
            "description": "Bases where overlapping mates of a pair disagree (tasmanian-diagnostics)",
            "scale": "Oranges",
            "format": "{:,.0f}",
            "hidden": True,
        }
    }
    headers = module.get_general_stats_headers(all_headers=headers)
    if headers:
        module.general_stats_addcols(
            {s_name: {"inconsistent_bases": d["inconsistent_bases"]} for s_name, d in inconsistency_data.items()},
            headers,
            namespace="tasmanian",
        )

    types = sorted({t for d in inconsistency_data.values() for t in d["by_type"]})
    add_category_bargraph(
        module,
        {s_name: {t: d["by_type"].get(t, 0) for t in types} for s_name, d in inconsistency_data.items()},
        name="Read Pair Inconsistencies",
        anchor="tasmanian-inconsistencies",
        description=(
            "Bases where the two mates of an overlapping read pair disagree, by the base seen in each mate, "
            "from <code>tasmanian-diagnostics</code>."
        ),
        helptext=(
            "Each category names the base observed in read 1 and in read 2 at the same genomic position, "
            "for example `R1:A_R2:G`."
        ),
        ylab="Bases",
    )

    module.write_data_file(inconsistency_data, "multiqc_tasmanian_inconsistencies")
    return len(inconsistency_data)


def parse_inconsistencies_table(contents: str) -> Dict:
    """Summarize a read_pair_inconsistencies.tsv table: read1_position, read2_position, discordance_type, count"""
    by_type: Dict[str, int] = {}
    for line in contents.splitlines()[1:]:
        if not line.strip():
            continue
        _read1_position, _read2_position, discordance_type, count = line.split("\t")
        if not DISCORDANCE_RE.match(discordance_type):
            raise ValueError(f"Unexpected discordance type: {discordance_type}")
        by_type[discordance_type] = by_type.get(discordance_type, 0) + int(count)
    return {"inconsistent_bases": sum(by_type.values()), "by_type": by_type}
