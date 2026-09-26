"""MultiQC submodule to parse potential_variants.tsv from tasmanian-diagnostics"""

import logging
from typing import Dict

from multiqc import BaseMultiqcModule

from .mismatch import MISMATCHES
from .util import add_category_bargraph

log = logging.getLogger(__name__)


def parse_tasmanian_variants(module: BaseMultiqcModule) -> int:
    """Find tasmanian-diagnostics potential variant tables and parse their data"""

    variants_data: Dict[str, Dict] = {}

    for f in module.find_log_files("tasmanian/variants"):
        try:
            data = parse_variants_table(f["f"])
        except ValueError as e:
            raise ValueError(f"Could not parse tasmanian-diagnostics output {f['fn']}: {e}") from e

        s_name = module.clean_s_name(f["s_name"], f)
        if s_name in variants_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        module.add_data_source(f, s_name=s_name, section="variants")
        variants_data[s_name] = data

        # tasmanian-diagnostics does not record its version in the output
        module.add_software_version(None, s_name)

    variants_data = module.ignore_samples(variants_data)
    if len(variants_data) == 0:
        return 0

    headers: Dict = {
        "variant_sites": {
            "title": "Potential Variant Sites",
            "description": "Genomic sites with recurrent mismatches (tasmanian-diagnostics)",
            "scale": "Purples",
            "format": "{:,.0f}",
            "hidden": True,
        }
    }
    headers = module.get_general_stats_headers(all_headers=headers)
    if headers:
        module.general_stats_addcols(
            {s_name: {"variant_sites": d["variant_sites"]} for s_name, d in variants_data.items()},
            headers,
            namespace="tasmanian",
        )

    add_category_bargraph(
        module,
        {s_name: d["sites_by_class"] for s_name, d in variants_data.items()},
        name="Potential Variant Sites",
        anchor="tasmanian-variants",
        description=(
            "Genomic sites where <code>tasmanian-diagnostics</code> found mismatches above its count and depth "
            "thresholds, by substitution class."
        ),
        helptext=(
            "Recurrent mismatches at one genomic position are more likely to be true variants than random artifacts."
        ),
        ylab="Sites",
    )

    module.write_data_file(variants_data, "multiqc_tasmanian_variants")
    return len(variants_data)


def parse_variants_table(contents: str) -> Dict:
    """Summarize a potential_variants.tsv table: chromosome, position, reference_base, mismatch_base, count, depth"""
    sites_by_class = {mismatch: 0 for mismatch in MISMATCHES}
    mismatch_count = 0
    for line in contents.splitlines()[1:]:
        if not line.strip():
            continue
        _chromosome, _position, ref_base, mismatch_base, count, _depth = line.split("\t")
        mismatch = f"{ref_base}>{mismatch_base}"
        if mismatch not in sites_by_class:
            raise ValueError(f"Unexpected substitution: {mismatch}")
        sites_by_class[mismatch] += 1
        mismatch_count += int(count)
    return {
        "variant_sites": sum(sites_by_class.values()),
        "mismatch_count": mismatch_count,
        "sites_by_class": sites_by_class,
    }
