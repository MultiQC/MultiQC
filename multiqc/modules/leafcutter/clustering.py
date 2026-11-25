import gzip
import logging
import os
import re
from operator import add, or_
from typing import Any, TextIO

from multiqc import config
from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph

log = logging.getLogger(__name__)

# Intron IDs are chrom:start:end:clu_<n>, with a _<strand> suffix from leafcutter_cluster_regtools.py
INTRON_ID_RE = re.compile(r"^(?P<chrom>\S+):\d+:\d+:(?P<cluster>clu_\d+(?:_\S+)?)$")


def parse_leafcutter_clustering(module: BaseMultiqcModule) -> int:
    """Parse the `*_perind_numers.counts.gz` junction count tables from LeafCutter clustering."""
    data: dict[str, dict[str, Any]] = {}

    for f in module.find_log_files("leafcutter/clustering", filecontents=False, filehandles=False):
        path = os.path.join(f["root"], f["fn"])
        # Matched on filename alone, as content search does not read gzip files
        try:
            with gzip.open(path, "rt") as fh:
                parsed = parse_perind_numers(fh)
        except (OSError, EOFError, UnicodeDecodeError, ValueError) as e:
            log.warning(f"Skipping {f['fn']}: {e}")
            continue

        for raw_name, sample_data in parsed.items():
            s_name = module.clean_s_name(raw_name, f)
            if s_name in data:
                log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {s_name}")
            data[s_name] = sample_data
            module.add_data_source(f, s_name=s_name, section="clustering")

    data = module.ignore_samples(data)
    if len(data) == 0:
        return 0

    module.add_software_version(None)

    module.general_stats_addcols(
        data,
        {
            "introns_detected_pct": {
                "title": "Introns Detected",
                "description": "Percentage of clustered introns with at least one supporting read",
                "min": 0,
                "max": 100,
                "suffix": "%",
                "scale": "RdYlGn",
            },
            "junction_reads": {
                "title": f"{config.read_count_prefix} Junction Reads",
                "description": (
                    f"Split reads on clustered introns ({config.read_count_desc}). A read spanning "
                    "several introns counts once for each, so this is not a read count."
                ),
                "min": 0,
                "scale": "Blues",
                "modify": lambda x: x * config.read_count_multiplier,
                "hidden": True,
            },
        },
    )

    module.add_section(
        name="Introns Detected",
        anchor="leafcutter_introns_detected",
        description="Clustered introns with and without a supporting read in each sample.",
        helptext="""
            LeafCutter pools the junctions of all samples, keeps introns with enough reads across
            the cohort and groups introns that share a splice site into clusters. Every sample is
            counted against the same set of introns, and an intron is detected in a sample when at
            least one read in that sample supports it.

            A sample that detects a much smaller share of the introns than the others has lower
            depth or library complexity. `leafcutter_ds` tests these counts directly, so such a
            sample contributes little and can push clusters below its coverage filters. The
            junction read counts behind this plot are in the General Statistics table (hidden by
            default).
            """,
        plot=bargraph.plot(
            data,
            {
                "introns_detected": {"name": "Detected (at least one read)", "color": "#437bb1"},
                "introns_not_detected": {"name": "Not detected", "color": "#b1b1b1"},
            },
            {
                "id": "leafcutter_introns_detected_plot",
                "title": "LeafCutter: Introns Detected",
                "ylab": "Introns",
                "cpswitch_counts_label": "Number of introns",
            },
        ),
    )

    module.write_data_file(data, "multiqc_leafcutter_clustering")

    return len(data)


def parse_perind_numers(fh: TextIO) -> dict[str, dict[str, int | float]]:
    """
    Parse a LeafCutter `*_perind_numers.counts.gz` table.

    The header lists the sample names, and each row starts with an intron ID followed by one
    read count per sample, all separated by single spaces.
    """
    header = fh.readline().split()
    if len(header) == 0:
        raise ValueError("no sample names in the header")

    junction_reads = [0] * len(header)
    introns_detected = [0] * len(header)
    # Per cluster, whether each sample has a read on any of its introns
    cluster_detected: dict[str, list[bool]] = {}
    n_introns = 0

    for line_number, line in enumerate(fh, start=2):
        fields = line.split()
        if not fields:
            continue
        if len(fields) != len(header) + 1:
            raise ValueError(f"expected {len(header) + 1} fields on line {line_number}, found {len(fields)}")
        match = INTRON_ID_RE.match(fields[0])
        if match is None:
            raise ValueError(f"unrecognised LeafCutter intron ID '{fields[0]}' on line {line_number}")
        cluster = f"{match.group('chrom')}:{match.group('cluster')}"
        counts = list(map(int, fields[1:]))
        detected = list(map(bool, counts))
        junction_reads = list(map(add, junction_reads, counts))
        introns_detected = list(map(add, introns_detected, detected))
        previous = cluster_detected.get(cluster)
        cluster_detected[cluster] = detected if previous is None else list(map(or_, previous, detected))
        n_introns += 1

    clusters_detected = [sum(column) for column in zip(*cluster_detected.values())] or [0] * len(header)

    parsed: dict[str, dict[str, int | float]] = {}
    for i, sample in enumerate(header):
        parsed[sample] = {
            "junction_reads": junction_reads[i],
            "introns": n_introns,
            "introns_detected": introns_detected[i],
            "introns_not_detected": n_introns - introns_detected[i],
            "clusters": len(cluster_detected),
            "clusters_detected": clusters_detected[i],
        }
        if n_introns > 0:
            parsed[sample]["introns_detected_pct"] = introns_detected[i] / n_introns * 100.0
    return parsed
