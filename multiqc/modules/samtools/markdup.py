import json
import logging
import re

from typing import Dict
from typing import Union

from multiqc.plots import bargraph
from multiqc.plots import table

log = logging.getLogger(__name__)


def parse_samtools_markdup(module) -> int:
    raw_by_sample: Dict = dict()

    for f in module.find_log_files("samtools/markdup_json", filehandles=True):
        raw_d = json.load(f["f"])
        if f["s_name"] in raw_by_sample:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {f['s_name']}")
        module.add_data_source(f, section="markdup")
        raw_by_sample[f["s_name"]] = raw_d

    for f in module.find_log_files("samtools/markdup_txt"):
        raw_d = dict()
        for line in f["f"].splitlines():
            if ":" in line:
                key, value = line.split(":")
                try:
                    value = int(value.strip())
                except ValueError:
                    value = value.strip()
                raw_d[key.strip()] = value
        if f["s_name"] in raw_by_sample:
            log.debug(f"Duplicate sample name found in {f['fn']}! Overwriting: {f['s_name']}")
        module.add_data_source(f, section="markdup")
        raw_by_sample[f["s_name"]] = raw_d

    raw_by_sample = module.ignore_samples(raw_by_sample)
    if len(raw_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    module.add_software_version(None)

    val_by_metric_by_sample: Dict[str, Dict[str, Union[int, float]]] = {}
    for s_name, raw_d in raw_by_sample.items():
        if len(raw_d) == 0:
            continue

        d: Dict[str, Union[int, float]] = {}
        val_by_metric_by_sample[s_name] = d
        rename_map = {
            "READ": "read",
            "WRITTEN": "written",
            "EXCLUDED": "excluded",
            "EXAMINED": "examined",
            "PAIRED": "paired",
            "SINGLE": "single",
            "DUPLICATE PAIR": "duplicate_pair",
            "DUPLICATE SINGLE": "duplicate_single",
            "DUPLICATE PAIR OPTICAL": "duplicate_pair_optical",
            "DUPLICATE SINGLE OPTICAL": "duplicate_single_optical",
            "DUPLICATE NON PRIMARY": "duplicate_non_primary",
            "DUPLICATE NON PRIMARY OPTICAL": "duplicate_non_primary_optical",
            "DUPLICATE PRIMARY TOTAL": "duplicate_primary_total",
            "DUPLICATE TOTAL": "duplicate_total",
            "ESTIMATED_LIBRARY_SIZE": "estimated_library_size",
        }
        for name, val in raw_d.items():
            if name in rename_map:
                d[rename_map[name]] = int(raw_d[name])
            elif name == "COMMAND":
                m = re.search(r".*-d\s(\d+)", val)
                if m:
                    d["optical_duplicate_distance"] = int(m.group(1))

        # Derive more metrics from counts
        n_reads: Union[int, float] = d["paired"] + d["single"]
        d["duplicate_optical_total"] = d["duplicate_pair_optical"] + d["duplicate_single_optical"]
        d["duplicate_optical_fraction"] = d["duplicate_optical_total"] / n_reads if n_reads > 0 else 0.0
        d["duplicate_fraction"] = d["duplicate_total"] / n_reads if n_reads > 0 else 0.0
        d["duplicate_pair_non_optical"] = d["duplicate_pair"] - d["duplicate_pair_optical"]
        d["duplicate_single_non_optical"] = d["duplicate_single"] - d["duplicate_single_optical"]
        d["duplicate_non_primary_non_optical"] = d["duplicate_non_primary"] - d["duplicate_non_primary_optical"]
        d["non_duplicate"] = d["paired"] + d["single"] - d["duplicate_total"]

    module.write_data_file(val_by_metric_by_sample, fn="multiqc_samtools_markdup")

    genstats_headers = {
        "duplicate_fraction": {
            "title": "Duplicates",
            "description": "The percent of all types of duplicate reads",
            "min": 0,
            "max": 100,
            "modify": lambda x: x * 100,
            "suffix": "%",
            "scale": "OrRd",
        },
        "estimated_library_size": {
            "title": "Est. library size",
            "description": "The estimated library size after de-duplication.",
            "min": 0,
            "format": "{:,d}",
        },
    }
    module.general_stats_addcols(data=val_by_metric_by_sample, headers=genstats_headers, namespace="markdup")

    module.add_section(
        name="Markdup: stats",
        anchor="samtools-markdup",
        description=(
            "Optical duplicates are due to either optical or clustering-based artifacts. "
            + "See the following links to learn more about instrument-based duplicate "
            + "artifacts:"
            + "<br>"
            + "<ul>"
            + '<li><a href="https://core-genomics.blogspot.com/2016/05/increased-read-duplication-on-patterned.html" '
            + 'target="_blank">Core Genomics Post: Increased Read Duplication on Patterned '
            + "Flowcells</a>"
            + "</li>"
            + '<li><a href="https://sequencing.qcfail.com/articles/illumina-patterned-flow-cells-generate-duplicated-sequences/" '
            + 'target="_blank">QC Fail Post: Illumina Patterned Flow Cells Generate '
            + "Duplicated Sequences</a>"
            + "</li>"
            + "</ul>."
        ),
        plot=table.plot(
            data=val_by_metric_by_sample,
            headers=dict(
                genstats_headers,
                **{
                    "optical_duplicate_distance": {
                        "title": "Optical distance",
                        "description": "The optical distance for considering instrument duplicates",
                        "min": 0,
                        "format": "{:,d}",
                        "scale": "RdYlGn",
                    },
                    "duplicate_optical_fraction": {
                        "title": "Optical dups",
                        "description": "The percent of optical/clustering duplicate reads",
                        "min": 0,
                        "max": 100,
                        "modify": lambda x: x * 100,
                        "suffix": "%",
                        "scale": "RdYlGn-rev",
                    },
                },
            ),
            pconfig={
                "id": "samtools-markdup-table",
                "title": "Samtools: duplicate-marked SAM records (alignments)",
            },
        ),
    )

    # Bar plot
    pconfig = {
        "id": "samtools-markdup-fraction",
        "title": "Samtools: markdup: duplicate categories",
        "ylab": "SAM Records",
    }

    keys: Dict[str, Dict[str, str]] = {
        "non_duplicate": {"name": "Non-duplicates"},
        "duplicate_pair_optical": {"name": "Optical duplicates in pairs"},
        "duplicate_single_optical": {"name": "Optical duplicates in singletons"},
        "duplicate_non_primary_optical": {"name": "Optical non-primary duplicate"},
        "duplicate_pair_non_optical": {"name": "Non-optical duplicates in pairs"},
        "duplicate_single_non_optical": {"name": "Non-optical duplicates in singletons"},
        "duplicate_non_primary_non_optical": {"name": "Non-optical non-primary duplicates"},
        "excluded": {"name": "Ignored (QC fail or unmapped)"},
    }

    module.add_section(
        name="Markdup: duplicate categories",
        anchor="samtools-markdup-categories",
        description=(
            "For more information about the duplicate categories, see the "
            + '<a href="https://www.htslib.org/doc/samtools-markdup.html#STATISTICS" '
            + 'target="_blank">samtools documentation</a>. '
        ),
        plot=bargraph.plot(data=val_by_metric_by_sample, cats=keys, pconfig=pconfig),
    )

    return len(val_by_metric_by_sample)
