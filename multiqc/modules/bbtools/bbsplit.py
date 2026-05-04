"""MultiQC submodule to parse BBSplit alignment statistics"""

import logging
from typing import Dict

from multiqc.base_module import BaseMultiqcModule
from multiqc.plots import bargraph, table

log = logging.getLogger(__name__)

# BBSplit stats file column order (after the reference name)
BBSPLIT_COLS = [
    "pct_unambiguous",
    "unambiguous_mb",
    "pct_ambiguous",
    "ambiguous_mb",
    "unambiguous_reads",
    "ambiguous_reads",
    "assigned_reads",
    "assigned_bases",
]


def parse_bbtools_bbsplit(module: BaseMultiqcModule) -> int:
    """
    Parse BBSplit output files showing read distribution across reference genomes.

    BBSplit aligns reads to multiple references simultaneously and determines
    which reference each read best matches.
    """

    data_by_sample: Dict[str, Dict] = {}

    for f in module.find_log_files("bbtools/bbsplit", filehandles=True):
        parsed = _parse_bbsplit_file(f["f"])
        if parsed:
            s_name = f["s_name"]
            if s_name in data_by_sample:
                log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
            module.add_data_source(f, s_name=s_name, section="bbsplit")
            data_by_sample[s_name] = parsed

    # Filter to strip out ignored sample names
    data_by_sample = module.ignore_samples(data_by_sample)

    if len(data_by_sample) == 0:
        return 0

    # Superfluous function call to confirm that it is used in this module
    module.add_software_version(None)

    # Collect all reference names across all samples
    all_ref_names = set()
    for sample_data in data_by_sample.values():
        all_ref_names.update(sample_data.keys())

    # Prepare data for general stats (assigned reads per reference)
    general_stats_data: Dict[str, Dict] = {}
    for s_name, sample_data in data_by_sample.items():
        general_stats_data[s_name] = {}
        for ref_name, ref_data in sample_data.items():
            general_stats_data[s_name][f"{ref_name}_assigned"] = ref_data.get("assigned_reads", 0)

    # General Stats headers
    general_stats_headers: Dict = {}
    for ref_name in sorted(all_ref_names):
        general_stats_headers[f"{ref_name}_assigned"] = {
            "title": f"{ref_name} Assigned",
            "description": f"Total number of reads assigned to {ref_name}",
            "scale": "Purples",
            "shared_key": "read_count",
        }

    general_stats_headers = module.get_general_stats_headers(all_headers=general_stats_headers)
    if general_stats_headers:
        module.general_stats_addcols(general_stats_data, general_stats_headers, namespace="bbsplit")

    # Prepare data for bar graph (assigned reads per reference)
    bargraph_data: Dict[str, Dict] = {}
    for s_name, sample_data in data_by_sample.items():
        bargraph_data[s_name] = {}
        for ref_name, ref_data in sample_data.items():
            bargraph_data[s_name][ref_name] = ref_data.get("assigned_reads", 0)

    # Create category definitions for bar graph
    cats = {ref_name: {"name": ref_name} for ref_name in sorted(all_ref_names)}

    pconfig = {
        "id": "bbtools-bbsplit-plot",
        "title": "BBTools: BBSplit Read Distribution",
        "ylab": "Number of Reads",
        "cpswitch_counts_label": "Number of Reads",
        "cpswitch_percent_label": "Percentage of Reads",
    }

    module.add_section(
        name="BBSplit Alignment Distribution",
        anchor="bbtools-bbsplit",
        description="Statistics on how many reads mapped to which reference genome.",
        helptext="Shows the percentage and count of reads aligned to each reference genome.",
        plot=bargraph.plot(bargraph_data, cats, pconfig),
    )

    # Create detailed table with all BBSplit columns
    table_data: Dict[str, Dict] = {}
    for s_name, sample_data in data_by_sample.items():
        table_data[s_name] = {}
        for ref_name, ref_data in sample_data.items():
            table_data[s_name][f"{ref_name}_assigned"] = ref_data.get("assigned_reads", 0)
            table_data[s_name][f"{ref_name}_unambig_pct"] = ref_data.get("pct_unambiguous", 0)
            table_data[s_name][f"{ref_name}_ambig_pct"] = ref_data.get("pct_ambiguous", 0)
            table_data[s_name][f"{ref_name}_unambig"] = ref_data.get("unambiguous_reads", 0)
            table_data[s_name][f"{ref_name}_ambig"] = ref_data.get("ambiguous_reads", 0)

    table_headers: Dict = {}
    for ref_name in sorted(all_ref_names):
        table_headers[f"{ref_name}_assigned"] = {
            "title": f"{ref_name} Assigned",
            "description": f"Total number of reads assigned to {ref_name}",
            "scale": "Purples",
            "shared_key": "read_count",
            "hidden": True,
        }
        table_headers[f"{ref_name}_unambig_pct"] = {
            "title": f"{ref_name} % Unambig",
            "description": f"Percentage of reads unambiguously aligned to {ref_name}",
            "suffix": "%",
            "scale": "Greens",
            "min": 0,
            "max": 100,
        }
        table_headers[f"{ref_name}_ambig_pct"] = {
            "title": f"{ref_name} % Ambig",
            "description": f"Percentage of reads ambiguously aligned to {ref_name}",
            "suffix": "%",
            "scale": "Oranges",
            "min": 0,
            "max": 100,
        }
        table_headers[f"{ref_name}_unambig"] = {
            "title": f"{ref_name} Unambig",
            "description": f"Number of reads unambiguously aligned to {ref_name}",
            "scale": "Blues",
            "shared_key": "read_count",
            "hidden": True,
        }
        table_headers[f"{ref_name}_ambig"] = {
            "title": f"{ref_name} Ambig",
            "description": f"Number of reads ambiguously aligned to {ref_name}",
            "scale": "Reds",
            "shared_key": "read_count",
            "hidden": True,
        }

    module.add_section(
        name="BBSplit Statistics",
        anchor="bbtools-bbsplit-stats",
        description="Detailed statistics on read alignment to each reference genome.",
        plot=table.plot(
            table_data,
            table_headers,
            pconfig={
                "id": "bbtools-bbsplit-stats-table",
                "namespace": "BBTools",
                "title": "BBTools: BBSplit Statistics",
            },
        ),
    )

    # Write parsed data to file
    module.write_data_file(data_by_sample, "multiqc_bbtools_bbsplit")

    return len(data_by_sample)


def _parse_bbsplit_file(f) -> Dict:
    """
    Parse a BBSplit output file.

    Expected columns: name, %unambiguousReads, unambiguousMB, %ambiguousReads,
    ambiguousMB, unambiguousReads, ambiguousReads, assignedReads, assignedBases
    """
    parsed_data: Dict = {}

    for line in f:
        line = line.strip()
        if not line:
            continue

        # Skip header line
        if line.startswith("#name"):
            continue

        parts = line.split("\t")
        if len(parts) < len(BBSPLIT_COLS) + 1:
            continue

        ref_name = parts[0]
        try:
            ref_data = {
                "pct_unambiguous": float(parts[1]),
                "unambiguous_mb": float(parts[2]),
                "pct_ambiguous": float(parts[3]),
                "ambiguous_mb": float(parts[4]),
                "unambiguous_reads": int(parts[5]),
                "ambiguous_reads": int(parts[6]),
                "assigned_reads": int(parts[7]),
                "assigned_bases": int(parts[8]),
            }
            parsed_data[ref_name] = ref_data
        except (ValueError, IndexError) as e:
            log.warning(f"Could not parse BBSplit line: {line} - {e}")
            continue

    return parsed_data
