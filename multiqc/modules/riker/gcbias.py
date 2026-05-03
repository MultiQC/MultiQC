"""Parse riker `gcbias` outputs (gcbias-detail.txt + gcbias-summary.txt)."""

import logging
from collections import defaultdict
from typing import Dict

from multiqc.plots import linegraph

from .util import read_tsv, to_float, to_int

log = logging.getLogger(__name__)

# Cap normalized coverage so the noisy high-GC tail doesn't compress the typical
# 0.5-1.5x range into a flat line. Matches riker's PDF chart (gcbias.rs).
_GCBIAS_Y_MAX = 2.0


def parse_reports(module):
    parsed_samples = set()
    parsed_samples |= _parse_summary(module)
    parsed_samples |= _parse_detail(module)
    return parsed_samples


def _parse_summary(module) -> set:
    data_by_sample: Dict[str, Dict[str, float]] = {}

    for f in module.find_log_files("riker/gcbias_summary", filehandles=True):
        for row in read_tsv(f["f"]):
            sample = row.pop("sample", None)
            if not sample:
                continue
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="gcbias_summary")
            data_by_sample[s_name] = {col: to_float(val) for col, val in row.items()}

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    # Dropout: lower is better. <2 is excellent, >5 concerning. Cap gradient at 10.
    headers = {
        "at_dropout": {
            "title": "AT dropout",
            "description": "Coverage deficit at GC 0-50% (lower is better)",
            "min": 0,
            "max": 10,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
        "gc_dropout": {
            "title": "GC dropout",
            "description": "Coverage deficit at GC 51-100% (lower is better)",
            "min": 0,
            "max": 10,
            "format": "{:,.2f}",
            "scale": "OrRd",
            "hidden": True,
        },
    }
    module.general_stats_addcols(data_by_sample, headers, namespace="gcbias")
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_gcbias_summary")
    return set(data_by_sample.keys())


def _parse_detail(module) -> set:
    # data_by_sample[sample][gc_pct] = normalized_coverage
    data_by_sample: Dict[str, Dict[int, float]] = defaultdict(dict)

    for f in module.find_log_files("riker/gcbias_detail", filehandles=True):
        rows_by_sample: Dict[str, Dict[int, float]] = defaultdict(dict)
        for row in read_tsv(f["f"]):
            sample = row.get("sample")
            if not sample:
                continue
            try:
                gc = to_int(row["gc"])
            except (KeyError, ValueError):
                continue
            rows_by_sample[sample][gc] = min(to_float(row.get("normalized_coverage", "")), _GCBIAS_Y_MAX)

        for sample, by_gc in rows_by_sample.items():
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="gcbias_detail")
            data_by_sample[s_name] = by_gc

    data_by_sample = dict(module.ignore_samples(data_by_sample))
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    pconfig = {
        "id": f"{module.anchor}_gcbias_normalized_coverage",
        "title": f"{module.name}: GC bias",
        "ylab": "Normalized coverage",
        "xlab": "GC content (%)",
        "x_decimals": False,
        "ymin": 0,
        "ymax": _GCBIAS_Y_MAX,
        "tt_label": "<b>{point.x}% GC</b>: {point.y:.3f}x",
    }
    module.add_section(
        name="GC bias",
        anchor=f"{module.anchor}-gcbias-detail",
        description=(
            "Normalized read coverage as a function of GC content. A flat line at 1.0 indicates no GC bias. "
            f"Values above {_GCBIAS_Y_MAX:g}x are clamped to keep the typical range readable."
        ),
        plot=linegraph.plot(data_by_sample, pconfig),
    )
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_gcbias_detail")
    return set(data_by_sample.keys())
