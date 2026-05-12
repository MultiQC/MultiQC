"""Parse riker `basic` outputs (base distribution, mean quality, quality score distribution)."""

import logging
from collections import defaultdict
from typing import Dict

from multiqc.plots import linegraph

from .util import read_tsv, to_int

_BASES = ("a", "c", "g", "t", "n")

log = logging.getLogger(__name__)


def parse_reports(module):
    parsed_samples = set()
    parsed_samples |= _parse_base_distribution(module)
    parsed_samples |= _parse_mean_quality(module)
    parsed_samples |= _parse_quality_distribution(module)
    return parsed_samples


def _parse_base_distribution(module) -> set:
    # data_by_sample[sample] = {cycle: (frac_a, frac_c, frac_g, frac_t, frac_n)}
    # Split read1 / read2 into separate samples so the plot shows two lines per input
    # (matching the picard module's BaseDistributionByCycle behaviour).
    data_by_sample: Dict[str, Dict[int, tuple]] = {}

    for f in module.find_log_files("riker/basic_base_dist", filehandles=True):
        rows_by_sample: Dict[str, Dict[int, Dict[int, tuple]]] = defaultdict(lambda: defaultdict(dict))
        for row in read_tsv(f["f"], source=f["fn"]):
            sample = row.get("sample")
            if not sample:
                continue
            try:
                read_end = to_int(row["read_end"])
                cycle = to_int(row["cycle"])
                tup = tuple(float(row[f"frac_{b}"]) * 100.0 for b in _BASES)
            except (KeyError, ValueError) as e:
                log.warning(f"riker: skipping row in {f['fn']} for sample {sample}: {e}")
                continue
            rows_by_sample[sample][read_end][cycle] = tup

        for sample, by_read_end in rows_by_sample.items():
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="basic_base_distribution")
            if 2 in by_read_end:
                data_by_sample[f"{s_name}_R1"] = by_read_end.get(1, {})
                data_by_sample[f"{s_name}_R2"] = by_read_end[2]
            else:
                data_by_sample[s_name] = by_read_end.get(1, {})

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    pconfig = {
        "id": f"{module.anchor}_base_distribution_by_cycle",
        "title": f"{module.name}: Base distribution by cycle",
        "ylab": "%",
        "xlab": "Cycle",
        "x_decimals": False,
        "tt_label": "<b>cycle {point.x}</b>: {point.y:.2f} %",
        "ymax": 100,
        "ymin": 0,
        "data_labels": [
            {"name": "% A", "ylab": "% A"},
            {"name": "% C", "ylab": "% C"},
            {"name": "% G", "ylab": "% G"},
            {"name": "% T", "ylab": "% T"},
            {"name": "% N", "ylab": "% N"},
        ],
    }
    series_data = [{} for _ in range(5)]
    for s_name, by_cycle in data_by_sample.items():
        for i in range(5):
            series_data[i][s_name] = {cycle: vals[i] for cycle, vals in by_cycle.items()}

    module.add_section(
        name="Base distribution by cycle",
        anchor=f"{module.anchor}-basic-base-dist",
        description="Per-cycle base composition from riker's `basic` tool. Read 1 and read 2 are shown as separate samples.",
        helptext="""
            Fraction of each base call (A, C, G, T, N) at every sequencing cycle, equivalent to
            Picard `CollectBaseDistributionByCycle`. Read 1 and read 2 are emitted as separate
            `_R1` / `_R2` samples so cycle 1 of each read starts at x = 1.

            Useful for spotting cycle-specific issues such as adapter readthrough at the ends of
            reads, base composition imbalance from a non-random library, or a failing channel
            visible as a sudden N spike.
        """,
        plot=linegraph.plot(series_data, pconfig),
    )
    flat: Dict[str, Dict[str, float]] = {}
    for s, by_cycle in data_by_sample.items():
        flat[s] = {}
        for cycle, vals in by_cycle.items():
            for i, base in enumerate(_BASES):
                flat[s][f"cycle_{cycle}_frac_{base}"] = vals[i] / 100.0
    module.write_data_file(flat, f"multiqc_{module.anchor}_basic_base_distribution")
    return {s.removesuffix("_R1").removesuffix("_R2") for s in data_by_sample.keys()}


def _parse_mean_quality(module) -> set:
    # Riker emits mean quality with cycles 1..N concatenated across read 1 and read 2.
    # We plot it as-is; no read_end split available.
    data_by_sample: Dict[str, Dict[int, float]] = {}

    for f in module.find_log_files("riker/basic_mean_quality", filehandles=True):
        rows_by_sample: Dict[str, Dict[int, float]] = defaultdict(dict)
        for row in read_tsv(f["f"], source=f["fn"]):
            sample = row.get("sample")
            if not sample:
                continue
            try:
                cycle = to_int(row["cycle"])
                rows_by_sample[sample][cycle] = float(row["mean_quality"])
            except (KeyError, ValueError) as e:
                log.warning(f"riker: skipping row in {f['fn']} for sample {sample}: {e}")
                continue

        for sample, by_cycle in rows_by_sample.items():
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="basic_mean_quality")
            data_by_sample[s_name] = by_cycle

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    pconfig = {
        "id": f"{module.anchor}_mean_quality_by_cycle",
        "title": f"{module.name}: Mean base quality by cycle",
        "ylab": "Mean base quality",
        "xlab": "Cycle",
        "x_decimals": False,
        "ymin": 0,
        "tt_label": "<b>cycle {point.x}</b>: {point.y:.2f}",
    }
    module.add_section(
        name="Mean base quality by cycle",
        anchor=f"{module.anchor}-basic-mean-quality",
        description="Mean base quality at each sequencing cycle. Cycles 1..N concatenate read 1 then read 2.",
        helptext="""
            Mean base quality (Phred) at each sequencing cycle, equivalent to Picard
            `MeanQualityByCycle`. Riker emits a single series per sample with read 1 cycles
            followed by read 2 cycles, so a paired-end run shows two humps with a quality
            reset at the boundary.

            A gradual decline towards the end of each read is normal. Sharp drops, plateaus
            near the lower bound, or a much weaker read 2 may indicate sequencing chemistry
            problems.
        """,
        plot=linegraph.plot(data_by_sample, pconfig),
    )
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_basic_mean_quality_by_cycle")
    return set(data_by_sample.keys())


def _parse_quality_distribution(module) -> set:
    data_by_sample: Dict[str, Dict[int, int]] = {}

    for f in module.find_log_files("riker/basic_quality_dist", filehandles=True):
        rows_by_sample: Dict[str, Dict[int, int]] = defaultdict(dict)
        for row in read_tsv(f["f"], source=f["fn"]):
            sample = row.get("sample")
            if not sample:
                continue
            try:
                quality = to_int(row["quality"])
                count = to_int(row["count"])
            except (KeyError, ValueError) as e:
                log.warning(f"riker: skipping row in {f['fn']} for sample {sample}: {e}")
                continue
            rows_by_sample[sample][quality] = count

        for sample, by_quality in rows_by_sample.items():
            s_name = module.clean_s_name(sample, f)
            module.add_data_source(f, s_name, section="basic_quality_distribution")
            data_by_sample[s_name] = by_quality

    data_by_sample = module.ignore_samples(data_by_sample)
    if not data_by_sample:
        return set()

    module.add_software_version(None)

    pconfig = {
        "id": f"{module.anchor}_quality_score_distribution",
        "title": f"{module.name}: Quality score distribution",
        "ylab": "Number of bases",
        "xlab": "Base quality",
        "x_decimals": False,
        "ymin": 0,
        "tt_label": "<b>Q{point.x}</b>: {point.y:,.0f}",
    }
    module.add_section(
        name="Quality score distribution",
        anchor=f"{module.anchor}-basic-quality-dist",
        description="Distribution of base quality scores across all bases.",
        helptext="""
            Histogram of base quality (Phred) scores aggregated across every base call in the BAM,
            equivalent to Picard `QualityScoreDistribution`. A healthy Illumina run typically shows
            most mass at Q30 and above. A heavy left tail or a bimodal distribution suggests
            quality issues that warrant a closer look at the per-cycle plot.
        """,
        plot=linegraph.plot(data_by_sample, pconfig),
    )
    module.write_data_file(data_by_sample, f"multiqc_{module.anchor}_basic_quality_score_distribution")
    return set(data_by_sample.keys())
