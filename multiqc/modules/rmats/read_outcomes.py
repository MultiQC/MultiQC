import logging
import re

from multiqc import BaseMultiqcModule, config
from multiqc.plots import bargraph
from multiqc.types import SectionAlert

log = logging.getLogger(__name__)

# The order and names rMATS-turbo writes, from READ_ENUM_NAMES in rmatspipeline.pyx
OUTCOMES: dict[str, dict[str, str]] = {
    "USED": {"name": "Used", "color": "#5cb85c"},
    "NOT_PAIRED": {"name": "Not properly paired", "color": "#fdae61"},
    "NOT_NH_1": {"name": "Multi-mapped (NH != 1)", "color": "#f46d43"},
    "NOT_EXPECTED_CIGAR": {"name": "Unexpected CIGAR", "color": "#d73027"},
    "NOT_EXPECTED_READ_LENGTH": {"name": "Unexpected read length", "color": "#a50026"},
    "NOT_EXPECTED_STRAND": {"name": "Unexpected strand", "color": "#8c510a"},
    "EXON_NOT_MATCHED_TO_ANNOTATION": {"name": "Exon not in annotation", "color": "#c7c7c7"},
    "JUNCTION_NOT_MATCHED_TO_ANNOTATION": {"name": "Junction not in annotation", "color": "#8a8a8a"},
    "CLIPPED": {"name": "Soft clipped", "color": "#fee08b"},
}

COUNT_LINE = re.compile(r"^([A-Z][A-Z0-9_]*): (\d+)$")

# Percentage of reads discarded by a default-off rMATS option above which the section warns
DISCARD_WARN_PCT = 1.0


def parse_read_outcomes_file(lines) -> dict[str, dict[str, int]]:
    """
    Parse one rMATS-turbo `read_outcomes_by_bam.txt` file. Each BAM is a block of a
    BAM path line, one `NAME: count` line per outcome, then `TOTAL_FOR_BAM: count`.

    Returns {bam_path: {outcome: count, "TOTAL_FOR_BAM": count}}.
    """
    blocks: dict[str, dict[str, int]] = {}
    bam: str | None = None
    for line in lines:
        line = line.rstrip("\r\n")
        if not line.strip():
            continue
        m = COUNT_LINE.match(line)
        if m is None:
            if bam is not None:
                raise ValueError(f"BAM '{bam}' has no TOTAL_FOR_BAM line before '{line}'")
            bam = line.strip()
            blocks[bam] = {}
            continue
        if bam is None:
            raise ValueError(f"count line '{line}' before any BAM path")
        key, value = m.group(1), int(m.group(2))
        blocks[bam][key] = value
        if key == "TOTAL_FOR_BAM":
            bam = None

    if bam is not None:
        raise ValueError(f"BAM '{bam}' has no TOTAL_FOR_BAM line")

    for bam_path, counts in blocks.items():
        missing = [k for k in OUTCOMES if k not in counts]
        if missing:
            raise ValueError(f"BAM '{bam_path}' is missing read outcomes: {', '.join(missing)}")
        summed = sum(v for k, v in counts.items() if k != "TOTAL_FOR_BAM")
        if summed != counts["TOTAL_FOR_BAM"]:
            raise ValueError(
                f"BAM '{bam_path}' outcomes sum to {summed}, but TOTAL_FOR_BAM is {counts['TOTAL_FOR_BAM']}"
            )
    return blocks


def bam_sample_name(module: BaseMultiqcModule, bam_path: str, f) -> str:
    """
    Sample name from the BAM path. When rMATS aligns FASTQ input itself, every BAM is
    `<tmp>/<prefix>_bam<group>_<replicate>/Aligned.sortedByCoord.out.bam`, so the
    directory name is used instead.
    """
    parts = [p for p in re.split(r"[\\/]", bam_path) if p]
    if parts[-1] == "Aligned.sortedByCoord.out.bam" and len(parts) > 1:
        return module.clean_s_name(parts[-2], f)
    return module.clean_s_name(parts[-1], f)


def parse_rmats_read_outcomes(module: BaseMultiqcModule) -> int:
    data: dict = {}

    for f in module.find_log_files("rmats/read_outcomes", filehandles=True):
        try:
            parsed = parse_read_outcomes_file(f["f"])
        except ValueError as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue
        for bam_path, counts in parsed.items():
            s_name = bam_sample_name(module, bam_path, f)
            if s_name in data and data[s_name] != counts:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            sample: dict[str, int | float] = dict(counts)
            total = float(counts["TOTAL_FOR_BAM"])
            if total > 0:
                sample["used_pct"] = counts["USED"] / total * 100.0
                sample["read_length_pct"] = counts["NOT_EXPECTED_READ_LENGTH"] / total * 100.0
                sample["clipped_pct"] = counts["CLIPPED"] / total * 100.0
            data[s_name] = sample
            module.add_data_source(f, s_name=s_name, section="read_outcomes")
            module.add_software_version(None, s_name)

    data = module.ignore_samples(data)
    if len(data) == 0:
        return 0

    general_stats_headers = {
        "used_pct": {
            "title": "Reads Used",
            "description": "Percentage of reads in the BAM that rMATS used for splicing event counts",
            "min": 0,
            "max": 100,
            "suffix": "%",
            "scale": "RdYlGn",
        },
        "TOTAL_FOR_BAM": {
            "title": f"{config.read_count_prefix} Reads",
            "description": f"Reads (BAM records) assessed by rMATS ({config.read_count_desc})",
            "min": 0,
            "scale": "Blues",
            "modify": lambda x: x * config.read_count_multiplier,
            "shared_key": "read_count",
            "hidden": True,
        },
    }
    module.general_stats_addcols(data, general_stats_headers)

    alerts = []
    no_reads_used = sorted(s for s, d in data.items() if d["USED"] == 0)
    if no_reads_used:
        alerts.append(
            SectionAlert(
                message="rMATS used no reads from these BAMs, so they contribute nothing to the event counts.",
                level="danger",
                affected_samples=no_reads_used,
            )
        )
    for key, message in [
        (
            "read_length_pct",
            (
                "reads were discarded for having a length other than `--readLength`. This is expected after "
                "adapter or quality trimming; rerun rMATS with `--variable-read-length` to keep them."
            ),
        ),
        (
            "clipped_pct",
            (
                "reads were discarded for being soft clipped. Aligners such as STAR soft clip by default; "
                "rerun rMATS with `--allow-clipping` to keep them."
            ),
        ),
    ]:
        affected = sorted(s for s, d in data.items() if d.get(key, 0.0) > DISCARD_WARN_PCT)
        if affected:
            alerts.append(
                SectionAlert(
                    message=f"In these samples more than {DISCARD_WARN_PCT:g}% of {message}",
                    level="warning",
                    affected_samples=affected,
                )
            )

    module.add_section(
        name="Read Outcomes",
        anchor="rmats_read_outcomes",
        description="How rMATS used the reads in each BAM file when counting splicing events.",
        helptext="""
        rMATS reads every record in each BAM file and either uses it to count exon and junction
        reads, or discards it for the first reason that applies:

        * **Soft clipped**: the alignment is soft clipped and `--allow-clipping` was not set.
        * **Not properly paired**: paired-end data where the read is not in a proper pair.
        * **Multi-mapped (NH != 1)**: the read does not map uniquely.
        * **Unexpected strand**: for a stranded library, the read is on the wrong strand for `-t`/`--libType`.
        * **Unexpected CIGAR**: the alignment has operations rMATS cannot use, such as insertions or deletions.
        * **Unexpected read length**: the aligned length differs from `--readLength`. Trimmed reads have
          variable lengths and are all discarded unless rMATS is run with `--variable-read-length`, which
          can silently remove most of the data.
        * **Exon not in annotation** / **Junction not in annotation**: the read has a usable layout but does
          not match an annotated exon or splice junction.

        A low proportion of used reads lowers the counts behind every splicing event, and so the power
        to detect differential splicing.
        """,
        alerts=alerts,
        plot=bargraph.plot(
            data,
            OUTCOMES,
            {
                "id": "rmats_read_outcomes_plot",
                "title": "rMATS: Read Outcomes",
                "ylab": "Reads",
                "cpswitch_counts_label": "Number of Reads",
            },
        ),
    )

    module.write_data_file(data, "multiqc_rmats_read_outcomes")

    return len(data)
