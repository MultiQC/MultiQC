import logging
import math
from collections.abc import Iterable

from multiqc import BaseMultiqcModule
from multiqc.plots import bargraph, linegraph
from multiqc.plots.bargraph import BarPlotConfig
from multiqc.plots.linegraph import LinePlotConfig
from multiqc.plots.table_object import ColumnDict

from .events import ISOFORM, event_type

log = logging.getLogger(__name__)

N_BINS = 20
BIN_WIDTH = 1.0 / N_BINS
PSI_TOLERANCE = 1e-9
LEVELS = {"event": "Local events", "isoform": "Isoforms"}


def parse_psi(lines: Iterable[str]) -> tuple[str, dict[str, dict]]:
    """
    Parse a SUPPA `.psi` matrix from `psiPerEvent` or `psiPerIsoform`.

    Returns the level (`event` or `isoform`) and per-sample summaries: the number of
    events, the number with a PSI value, and a histogram of PSI values.
    """
    line_iter = iter(lines)
    try:
        samples = next(line_iter).rstrip("\r\n").split("\t")
    except StopIteration:
        raise ValueError("file is empty")
    if not all(samples):
        raise ValueError("header has empty sample names")
    if len(set(samples)) != len(samples):
        raise ValueError("header has duplicate sample names")

    quantified = dict.fromkeys(samples, 0)
    psi_sums = dict.fromkeys(samples, 0.0)
    bins = {s: [0] * N_BINS for s in samples}
    n_events = 0
    level = None
    for line_no, line in enumerate(line_iter, start=2):
        fields = line.rstrip("\r\n").split("\t")
        if fields == [""]:
            continue
        if len(fields) != len(samples) + 1:
            raise ValueError(f"line {line_no} has {len(fields) - 1} values for {len(samples)} samples")
        row_level = "isoform" if event_type(fields[0]) == ISOFORM else "event"
        if level is None:
            level = row_level
        elif row_level != level:
            raise ValueError(f"line {line_no} mixes local event and isoform ids")
        n_events += 1
        for sample, value in zip(samples, fields[1:]):
            # SUPPA writes `nan` or `NA` when PSI cannot be computed. Releases before 2.0
            # wrote -1, so any negative value is treated the same way.
            if value == "NA":
                continue
            psi = float(value)
            if math.isnan(psi) or psi < 0:
                continue
            # Floating point sums can give values such as 1.0000000000000002
            if psi > 1 + PSI_TOLERANCE:
                raise ValueError(f"line {line_no} has PSI {value} outside [0, 1]")
            quantified[sample] += 1
            psi_sums[sample] += psi
            bins[sample][min(int(psi / BIN_WIDTH), N_BINS - 1)] += 1

    if level is None:
        raise ValueError("file has no events")

    data: dict[str, dict] = {}
    for sample in samples:
        n = quantified[sample]
        data[sample] = {
            "events": n_events,
            "quantified": n,
            "not_quantified": n_events - n,
            "pct_quantified": 100.0 * n / n_events,
        }
        if n:
            data[sample]["mean_psi"] = psi_sums[sample] / n
            data[sample]["histogram"] = {
                round((i + 0.5) * BIN_WIDTH, 3): 100.0 * count / n for i, count in enumerate(bins[sample])
            }
    return level, data


def parse_suppa_psi(module: BaseMultiqcModule) -> int:
    """Find SUPPA `.psi` files and add PSI sections. Returns the number of samples."""
    data: dict[str, dict[str, dict]] = {level: {} for level in LEVELS}

    for f in module.find_log_files("suppa/psi", filehandles=True):
        try:
            level, parsed = parse_psi(f["f"])
        except ValueError as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue
        for sample, sample_data in parsed.items():
            s_name = module.clean_s_name(sample, f)
            if s_name in data[level]:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section=f"psi_{level}")
            module.add_software_version(None, s_name)
            data[level][s_name] = sample_data

    for level in LEVELS:
        data[level] = module.ignore_samples(data[level])
    levels = [level for level in LEVELS if data[level]]
    if not levels:
        return 0

    _add_general_stats(module, data)
    _add_quantified_section(module, data, levels)
    _add_distribution_section(module, data, levels)

    for level in levels:
        module.write_data_file(
            {s: {k: v for k, v in d.items() if k != "histogram"} for s, d in data[level].items()},
            f"multiqc_suppa_psi_{level}",
        )
    return len(set().union(*(data[level] for level in levels)))


def _add_general_stats(module: BaseMultiqcModule, data: dict[str, dict[str, dict]]) -> None:
    for level, label, hidden in [("event", "events", False), ("isoform", "isoforms", True)]:
        if not data[level]:
            continue
        headers: dict[str, ColumnDict] = {
            f"pct_quantified_{level}": {
                "title": f"{label.capitalize()} Quantified",
                "description": f"Percentage of SUPPA {label} with a PSI value (not NA or nan)",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
                "hidden": hidden,
            }
        }
        module.general_stats_addcols(
            {s: {f"pct_quantified_{level}": d["pct_quantified"]} for s, d in data[level].items()},
            headers,
        )


def _add_quantified_section(module: BaseMultiqcModule, data: dict[str, dict[str, dict]], levels: list[str]) -> None:
    cats = {
        "quantified": {"name": "Quantified (PSI value)", "color": "#437bb1"},
        "not_quantified": {"name": "Not quantified (NA or nan)", "color": "#b1b1b1"},
    }
    module.add_section(
        name="PSI Quantification",
        anchor="suppa_psi_quantification",
        description="Number of local splicing events (`psiPerEvent`) or isoforms (`psiPerIsoform`) "
        "for which SUPPA could compute a PSI value in each sample.",
        helptext="""
SUPPA reports `nan` or `NA` instead of a PSI value when the transcripts that define an
event are not expressed in a sample, when their total expression is below the
`--total-filter` threshold, or when a transcript is missing from the expression file.

A sample with a much lower fraction of quantified events than others usually has low
depth or low transcript-level expression, and will contribute fewer events to
differential splicing tests (by default `diffSplice` discards any event with a missing
value in any replicate, see `--nan-threshold`).
""",
        plot=bargraph.plot(
            [{s: {k: d[k] for k in cats} for s, d in data[level].items()} for level in levels],
            [cats] * len(levels),
            BarPlotConfig(
                id="suppa_psi_quantification_plot",
                title="SUPPA: PSI Quantification",
                ylab="Events",
                cpswitch_counts_label="Number of events",
                data_labels=[{"name": LEVELS[level], "ylab": LEVELS[level]} for level in levels],
            ),
        ),
    )


def _add_distribution_section(module: BaseMultiqcModule, data: dict[str, dict[str, dict]], levels: list[str]) -> None:
    module.add_section(
        name="PSI Distribution",
        anchor="suppa_psi_distribution",
        description=f"Distribution of PSI values in each sample, in bins of {BIN_WIDTH:g}, "
        "as a percentage of the events with a PSI value.",
        helptext="""
PSI (percent spliced in) is the fraction of the expression of an event's transcripts
that comes from the transcripts including the alternative form. Most local events and
isoforms sit close to 0 or 1, so the distribution is usually U-shaped. Samples from the
same condition should have similar curves; an outlier can point to a quantification
problem in that sample.
""",
        plot=linegraph.plot(
            [{s: d["histogram"] for s, d in data[level].items() if "histogram" in d} for level in levels],
            LinePlotConfig(
                id="suppa_psi_distribution_plot",
                title="SUPPA: PSI Distribution",
                xlab="PSI",
                ylab="% of quantified events",
                xmin=0,
                xmax=1,
                ymin=0,
                tt_suffix="%",
                data_labels=[
                    {"name": LEVELS[level], "ylab": f"% of quantified {LEVELS[level].lower()}"} for level in levels
                ],
            ),
        ),
    )
