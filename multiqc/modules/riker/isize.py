"""Parse riker `isize` outputs (isize-metrics.txt + isize-histogram.txt)."""

import logging
from collections import defaultdict
from typing import Dict

from multiqc.plots import linegraph

from .util import read_tsv, to_float, to_int

log = logging.getLogger(__name__)


def parse_reports(module):
    metrics_by_sample = _parse_metrics(module)
    parsed_samples = set(metrics_by_sample.keys())
    parsed_samples |= _parse_histogram(module)
    return parsed_samples


def _parse_metrics(module) -> Dict[str, Dict[str, Dict[str, float]]]:
    # data_by_sample[sample][orientation] = {col: value}
    data_by_sample: Dict[str, Dict[str, Dict[str, float]]] = {}

    for f in module.find_log_files("riker/isize_metrics", filehandles=True):
        for row in read_tsv(f["f"]):
            sample = row.pop("sample", None)
            orientation = row.pop("pair_orientation", None)
            if not sample or not orientation:
                continue
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="isize_metrics")
            parsed = {col: to_float(val) for col, val in row.items()}
            data_by_sample.setdefault(s_name, {})[orientation] = parsed

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return {}

    module.add_software_version(None)

    # Pick the FR orientation when present (the typical Illumina library), else the
    # first available orientation — there is always at least one row per sample.
    primary_data: Dict[str, Dict[str, float]] = {}
    for s_name, by_orient in data_by_sample.items():
        primary_data[s_name] = by_orient.get("FR", next(iter(by_orient.values())))

    # Insert size is library-design dependent: there is no universally "good"
    # value. Use a neutral single-hue scale so cells convey magnitude only.
    headers = {
        "median_insert_size": {
            "title": "Median insert",
            "description": "Median insert size in bp (FR orientation when present)",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
        },
        "mean_insert_size": {
            "title": "Mean insert",
            "description": "Mean insert size in bp (FR orientation when present)",
            "min": 0,
            "suffix": " bp",
            "format": "{:,.0f}",
            "scale": "GnBu",
            "hidden": True,
        },
    }
    module.general_stats_addcols(primary_data, headers, namespace="isize")
    module.write_data_file(primary_data, f"multiqc_{module.anchor}_isize")
    return data_by_sample


def _parse_histogram(module) -> set:
    # Per-orientation histogram. Keep the three orientations as separate datasets.
    fr_by_sample: Dict[str, Dict[int, int]] = defaultdict(dict)
    rf_by_sample: Dict[str, Dict[int, int]] = defaultdict(dict)
    tandem_by_sample: Dict[str, Dict[int, int]] = defaultdict(dict)

    for f in module.find_log_files("riker/isize_histogram", filehandles=True):
        rows_by_sample: Dict[str, list] = defaultdict(list)
        for row in read_tsv(f["f"]):
            sample = row.get("sample")
            if not sample:
                continue
            try:
                isize = to_int(row["insert_size"])
                fr = to_int(row.get("fr_count", "0"))
                rf = to_int(row.get("rf_count", "0"))
                tandem = to_int(row.get("tandem_count", "0"))
            except (KeyError, ValueError):
                continue
            rows_by_sample[sample].append((isize, fr, rf, tandem))

        for sample, rows in rows_by_sample.items():
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="isize_histogram")
            for isize, fr, rf, tandem in rows:
                if fr:
                    fr_by_sample[s_name][isize] = fr
                if rf:
                    rf_by_sample[s_name][isize] = rf
                if tandem:
                    tandem_by_sample[s_name][isize] = tandem

    fr_by_sample = dict(module.ignore_samples(fr_by_sample))
    rf_by_sample = dict(module.ignore_samples(rf_by_sample))
    tandem_by_sample = dict(module.ignore_samples(tandem_by_sample))

    sample_set = set(fr_by_sample.keys()) | set(rf_by_sample.keys()) | set(tandem_by_sample.keys())
    if not sample_set:
        return set()

    module.add_software_version(None)

    pconfig = {
        "id": f"{module.anchor}_insert_size",
        "title": f"{module.name}: Insert size distribution",
        "ylab": "Read pair count",
        "xlab": "Insert size (bp)",
        "x_decimals": False,
        "ymin": 0,
        "smooth_points": 500,
        "smooth_points_sumcounts": True,
        "tt_label": "<b>{point.x} bp</b>: {point.y:,.0f}",
        "data_labels": [
            {"name": "FR", "ylab": "Read pair count"},
            {"name": "RF", "ylab": "Read pair count"},
            {"name": "Tandem", "ylab": "Read pair count"},
        ],
    }
    module.add_section(
        name="Insert size distribution",
        anchor=f"{module.anchor}-isize-histogram",
        description="Insert size distribution per read-pair orientation.",
        plot=linegraph.plot([fr_by_sample, rf_by_sample, tandem_by_sample], pconfig),
    )
    return sample_set
