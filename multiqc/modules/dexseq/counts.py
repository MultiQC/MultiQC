import logging
from collections.abc import Iterable

from multiqc import BaseMultiqcModule, config
from multiqc.plots import bargraph
from multiqc.plots.table_object import ColumnDict

log = logging.getLogger(__name__)

# Counters that every version of dexseq_count.py writes, mapped to their data keys
REQUIRED_COUNTERS = {
    "_empty": "empty",
    "_ambiguous": "ambiguous",
    "_lowaqual": "lowaqual",
    "_notaligned": "notaligned",
}
# Absent from dexseq_count.py output written before 2014
OPTIONAL_COUNTERS = {
    "_ambiguous_readpair_position": "ambiguous_readpair_position",
}


def parse_counts_file(lines: Iterable[str]) -> dict[str, int | float] | None:
    """
    Summarise one dexseq_count.py output file.

    Returns None when the file holds none of the special counters, meaning it is not
    dexseq_count.py output. Raises ValueError when it holds some of them but not all,
    or has a malformed line, as that means the format has changed.
    """
    counters: dict[str, int] = {}
    bins_total = 0
    bins_detected = 0
    bin_counts = 0
    for line_number, line in enumerate(lines, start=1):
        line = line.rstrip("\r\n")
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) != 2:
            raise ValueError(f"line {line_number}: expected 2 tab-separated fields, found {len(fields)}")
        name, value = fields
        try:
            count = int(value)
        except ValueError:
            raise ValueError(f"line {line_number}: count '{value}' is not an integer")
        if name in REQUIRED_COUNTERS:
            counters[REQUIRED_COUNTERS[name]] = count
        elif name in OPTIONAL_COUNTERS:
            counters[OPTIONAL_COUNTERS[name]] = count
        else:
            bins_total += 1
            bin_counts += count
            if count > 0:
                bins_detected += 1

    if not counters:
        return None
    missing = [name for name, key in REQUIRED_COUNTERS.items() if key not in counters]
    if missing:
        raise ValueError(f"missing dexseq_count.py counters: {', '.join(missing)}")

    data: dict[str, int | float] = dict(counters)
    data["bins_total"] = bins_total
    data["bins_detected"] = bins_detected
    data["bins_not_detected"] = bins_total - bins_detected
    data["bin_counts"] = bin_counts
    data["unassigned"] = sum(counters.values())
    if bins_total > 0:
        data["bins_detected_pct"] = bins_detected / bins_total * 100.0
    return data


def parse_dexseq_counts(module: BaseMultiqcModule) -> int:
    """Find dexseq_count.py outputs, add their sections and return the number of samples"""
    data_by_sample: dict[str, dict[str, int | float]] = {}
    for f in module.find_log_files("dexseq/counts", filehandles=True):
        try:
            parsed = parse_counts_file(f["f"])
        except ValueError as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue
        if parsed is None:
            log.debug(f"Skipping {f['fn']}: no dexseq_count.py counters found")
            continue
        s_name = f["s_name"]
        if s_name in data_by_sample:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        module.add_data_source(f, s_name=s_name, section="counts")
        module.add_software_version(None, s_name)
        data_by_sample[s_name] = parsed

    data_by_sample = module.ignore_samples(data_by_sample)
    if len(data_by_sample) == 0:
        return 0

    _add_general_stats(module, data_by_sample)
    _add_bins_detected_section(module, data_by_sample)

    module.write_data_file(data_by_sample, "multiqc_dexseq_counts")
    return len(data_by_sample)


def _add_general_stats(module: BaseMultiqcModule, data_by_sample: dict[str, dict[str, int | float]]) -> None:
    headers: dict[str, ColumnDict] = {
        "bins_detected_pct": {
            "title": "Bins Detected",
            "description": "Percentage of exon bins with at least one supporting read",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
        },
        "unassigned": {
            "title": f"{config.read_count_prefix} Unassigned",
            "description": (
                f"Reads or read pairs not counted to any exon bin ({config.read_count_desc}). Compare with "
                "the aligned reads: a large share with no feature often means the wrong --stranded setting."
            ),
            "min": 0,
            "scale": "OrRd",
            "modify": lambda x: float(x) * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        },
    }
    module.general_stats_addcols({s_name: dict(data) for s_name, data in data_by_sample.items()}, headers)


def _add_bins_detected_section(module: BaseMultiqcModule, data_by_sample: dict[str, dict[str, int | float]]) -> None:
    cats = {
        "bins_detected": {"name": "Detected (at least one read)", "color": "#437bb1"},
        "bins_not_detected": {"name": "Not detected", "color": "#b1b1b1"},
    }
    module.add_section(
        name="Exon Bin Detection",
        anchor="dexseq_exon_bin_detection",
        description="Exon bins of the flattened annotation with and without a supporting read in each sample.",
        helptext="""
            `dexseq_prepare_annotation.py` splits the exons of each gene into non-overlapping
            exon bins (exonic parts). Every sample is counted against the same bins, and a bin is
            detected when it has at least one count. A sample that detects a much smaller share of
            the bins than the others usually has lower depth or a library problem.

            The unassigned read counters that `dexseq_count.py` writes are in the General
            Statistics table (hidden by default) and the module data file. They are not plotted:
            a read overlapping several exon bins of one gene is counted once for each bin, and the
            script does not write how many reads it counted, so there is no total to put them
            against. Compare them with the aligned reads of the aligner instead.
            """,
        plot=bargraph.plot(
            data_by_sample,
            cats,
            {
                "id": "dexseq_exon_bin_detection_plot",
                "title": "DEXSeq: Exon Bin Detection",
                "ylab": "Exon bins",
                "cpswitch_counts_label": "Number of exon bins",
            },
        ),
    )
