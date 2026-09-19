"""MultiQC submodule to parse output from tasmanian-mismatch"""

import logging
from collections import defaultdict
from typing import Dict, List, Mapping, Optional

from multiqc import BaseMultiqcModule, config
from multiqc.plots import heatmap, linegraph, table
from multiqc.utils import mqc_colour

log = logging.getLogger(__name__)

BASES = "ACGT"
MISMATCHES = [f"{ref}>{alt}" for ref in BASES for alt in BASES if ref != alt]
ALL_BASE_CHANGES = {f"{ref}>{alt}" for ref in BASES for alt in BASES}
ALL_MISMATCHES = "All Mismatches"

# Pooled classes: the rate is the summed count of the member substitutions over all bases with
# one of the member reference bases, so a group is not an average of its members' rates
GROUPS: Dict[str, List[str]] = {
    "Deamination (C>T + G>A)": ["C>T", "G>A"],
    "Oxidation (G>T + C>A)": ["G>T", "C>A"],
}
# Class name -> the substitutions it counts, for the profile plot, heatmaps and summary table
CLASSES: Dict[str, List[str]] = {ALL_MISMATCHES: MISMATCHES, **GROUPS, **{m: [m] for m in MISMATCHES}}
HEATMAP_CLASSES = [ALL_MISMATCHES, *GROUPS]
HEATMAP_UNITS = "Mismatches (%)"

# read number -> position -> base change -> count
Profile = Dict[int, Dict[int, Dict[str, int]]]


def parse_tasmanian_mismatch(module: BaseMultiqcModule) -> int:
    """Find tasmanian-mismatch tables and parse their data"""

    mismatch_data: Dict[str, Dict] = {}

    for f in module.find_log_files("tasmanian/mismatch"):
        try:
            data = parse_mismatch_table(f["f"])
        except (ValueError, IndexError) as e:
            raise ValueError(f"Could not parse tasmanian-mismatch output {f['fn']}: {e}") from e

        s_name = module.clean_s_name(f["s_name"], f)
        if s_name in mismatch_data:
            log.debug(f"Duplicate sample name found! Overwriting: {s_name}")
        module.add_data_source(f, s_name=s_name, section="mismatch")
        mismatch_data[s_name] = data

        # tasmanian-mismatch does not record its version in the output
        module.add_software_version(None, s_name)

    mismatch_data = module.ignore_samples(mismatch_data)
    if len(mismatch_data) == 0:
        return 0

    _add_general_stats(module, mismatch_data)
    _add_summary_table(module, mismatch_data)
    _add_profile_section(module, mismatch_data)
    _add_heatmap_sections(module, mismatch_data)

    module.write_data_file(mismatch_data, "multiqc_tasmanian_mismatch")
    return len(mismatch_data)


def parse_mismatch_table(contents: str) -> Dict:
    """
    Parse a tasmanian-mismatch count table (insert position mode) into a profile keyed by read number
    and fragment position. Counts are summed over `reference_order`.
    """
    lines = contents.splitlines()
    header = lines[0].split("\t")
    if header[3] != "fragment_position":
        raise ValueError(f"Unexpected position column: {header[3]}")
    if header[4] != "count":
        raise ValueError(f"Unexpected value column: {header[4]}")

    profile: Profile = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    for line in lines[1:]:
        if not line.strip():
            continue
        base_change, read_num, _reference_order, position, value = line.split("\t")
        if base_change not in ALL_BASE_CHANGES:
            raise ValueError(f"Unexpected base change: {base_change}")
        profile[int(read_num)][int(position)][base_change] += int(value)

    return {
        "profile": {
            rn: {pos: dict(bcs) for pos, bcs in sorted(positions.items())} for rn, positions in profile.items()
        },
    }


def class_rate(counts: Mapping[str, float], members: List[str]) -> Optional[float]:
    """
    Percentage of bases with a reference base of any member that carry one of the member substitutions.
    None if there are no such bases.
    """
    ref_bases = {m[0] for m in members}
    total = sum(counts.get(f"{ref}>{alt}", 0) for ref in ref_bases for alt in BASES)
    if total == 0:
        return None
    return 100.0 * sum(counts.get(m, 0) for m in members) / total


def overall_rates(data: Dict) -> Optional[Dict[str, float]]:
    """
    Overall mismatch rate per class (see `CLASSES`), over all positions and reads, plus the number of bases
    counted as "total_bases". Returns None for a table without any counts.
    """
    totals: Dict[str, float] = defaultdict(float)
    for positions in data["profile"].values():
        for base_changes in positions.values():
            for base_change, count in base_changes.items():
                totals[base_change] += count

    total_bases = sum(totals.values())
    if total_bases == 0:
        return None
    rates: Dict[str, float] = {"total_bases": total_bases}
    for name, members in CLASSES.items():
        rate = class_rate(totals, members)
        if rate is not None:
            rates[name] = rate
    return rates


def _add_general_stats(module: BaseMultiqcModule, mismatch_data: Dict[str, Dict]) -> None:
    stats: Dict = {}
    for s_name, data in mismatch_data.items():
        rates = overall_rates(data)
        if rates is not None:
            stats[s_name] = {"mismatch_rate": rates[ALL_MISMATCHES], "total_bases": rates["total_bases"]}
    headers: Dict = {
        "mismatch_rate": {
            "title": "Mismatch Rate",
            "description": "Percentage of aligned bases that differ from the reference (tasmanian-mismatch)",
            "min": 0,
            "suffix": "%",
            "scale": "OrRd",
            "format": "{:,.3f}",
        },
        "total_bases": {
            "title": "Bases Analyzed",
            "description": f"Aligned bases counted by tasmanian-mismatch ({config.base_count_desc})",
            "scale": "Greens",
            "shared_key": "base_count",
            "hidden": True,
        },
    }
    headers = module.get_general_stats_headers(all_headers=headers)
    if stats and headers:
        module.general_stats_addcols(stats, headers, namespace="tasmanian")


def _add_summary_table(module: BaseMultiqcModule, mismatch_data: Dict[str, Dict]) -> None:
    rates = {s_name: r for s_name, data in mismatch_data.items() if (r := overall_rates(data)) is not None}
    if not rates:
        return
    headers: Dict = {}
    for name, members in CLASSES.items():
        if name == ALL_MISMATCHES:
            description = "Percentage of all aligned bases that mismatch"
        else:
            refs = " or ".join(sorted({m[0] for m in members}))
            description = f"Percentage of bases with reference base {refs} that mismatch as {' or '.join(members)}"
        headers[name] = {
            "title": name,
            "description": description,
            "min": 0,
            "suffix": "%",
            "scale": "OrRd",
            "format": "{:,.3f}",
            "hidden": name in MISMATCHES,
        }

    module.add_section(
        name="Mismatch Rates",
        anchor="tasmanian-mismatch-rates",
        description=(
            "Overall mismatch rates from <code>tasmanian-mismatch</code>, over all positions. "
            "Each substitution class is the percentage of bases with that reference base "
            "that were sequenced as the alternative base. Deamination and oxidation pool their two "
            "substitutions over all bases with either reference base. Individual substitutions are "
            "hidden by default; show them with the columns button."
        ),
        plot=table.plot(
            rates,
            headers,
            pconfig={
                "id": "tasmanian-mismatch-rates-table",
                "title": "Tasmanian: Mismatch Rates",
                "namespace": "tasmanian",
            },
        ),
    )


def _count_rates(positions: Dict[int, Dict[str, float]]) -> Dict[str, Dict[int, float]]:
    """Percentages by class (see `CLASSES`) and position, from raw counts."""
    by_class: Dict[str, Dict[int, float]] = {name: {} for name in CLASSES}
    for pos, counts in positions.items():
        for name, members in CLASSES.items():
            rate = class_rate(counts, members)
            if rate is not None:
                by_class[name][pos] = rate
    return by_class


def _profile_rates(data: Dict) -> Dict[str, Dict[int, float]]:
    """Convert a parsed profile to percentages by class and position, with the reads merged"""
    merged: Dict[int, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
    for positions in data["profile"].values():
        for pos, counts in positions.items():
            for base_change, count in counts.items():
                merged[pos][base_change] += count
    return _count_rates(merged)


def _add_profile_section(module: BaseMultiqcModule, mismatch_data: Dict[str, Dict]) -> None:
    datasets: List[Dict[str, Dict[int, float]]] = [{} for _ in CLASSES]
    # Color by sample, so that a trace keeps its color when another substitution class is chosen
    scale = mqc_colour.mqc_colour_scale("plot_defaults")
    colors: Dict[str, str] = {}
    for sample_idx, (s_name, data) in enumerate(mismatch_data.items()):
        colors[s_name] = scale.get_colour(sample_idx, lighten=1)
        by_class = _profile_rates(data)
        for i, name in enumerate(CLASSES):
            if by_class[name]:
                datasets[i][s_name] = by_class[name]

    data_labels = [{"name": name, "ylab": f"{name} (%)"} for name in CLASSES]
    # Drop classes without any bases, such as a reference base that never occurs
    keep = [i for i, ds in enumerate(datasets) if ds]
    datasets, data_labels = [datasets[i] for i in keep], [data_labels[i] for i in keep]

    module.add_section(
        name="Mismatch Profile by Normalized Fragment Position",
        anchor="tasmanian-profile",
        description=(
            "Percentage of bases that mismatch the reference at each normalized fragment position, "
            "from <code>tasmanian-mismatch</code>. Choose a substitution class with the buttons above the plot; "
            "each class is a percentage of the bases with that reference base."
        ),
        helptext=(
            "Artifacts such as oxidation (G>T and C>A) or deamination (C>T and G>A) appear as raised "
            "rates concentrated at the ends of fragments, whereas true variants and sequencing "
            "errors are spread evenly along the fragment. The deamination and oxidation buttons pool "
            "their two substitutions, and read 1 and read 2 are merged."
        ),
        plot=linegraph.plot(
            datasets,
            pconfig={
                "id": "tasmanian-profile-plot",
                "title": "Tasmanian: Mismatch Profile by Normalized Fragment Position",
                "xlab": "Normalized Fragment Position",
                "ylab": "Mismatches (%)",
                "ymin": 0,
                "tt_suffix": "%",
                "tt_decimals": 3,
                "colors": colors,
                "data_labels": data_labels,
            },
        ),
    )


def _add_heatmap_sections(module: BaseMultiqcModule, mismatch_data: Dict[str, Dict]) -> None:
    """One heatmap per pooled class, with a row per library, to compare many libraries at once"""
    rates = {s_name: _profile_rates(d) for s_name, d in mismatch_data.items()}

    for name in HEATMAP_CLASSES:
        data: Dict = {s_name: by_class[name] for s_name, by_class in rates.items()}
        positions = sorted({pos for by_pos in data.values() for pos in by_pos})
        slug = name.split(" (")[0].lower().replace(" ", "-")
        plot = heatmap.plot(
            data,
            xcats=positions,
            ycats=list(data),
            pconfig={
                "id": f"tasmanian-heatmap-{slug}-plot",
                "title": f"Tasmanian: {name}",
                "xlab": "Normalized Fragment Position",
                "ylab": "",
                "zlab": HEATMAP_UNITS,
                "min": 0,
                "tt_decimals": 3,
                "xcats_samples": False,
                "square": False,
                "cluster_cols": False,
                "angled_xticks": False,
            },
        )
        if isinstance(plot, heatmap.HeatmapPlot):
            # The heatmap module does not label the color scale, so add the units to the legend
            for dataset in plot.datasets:
                dataset.trace_params["colorbar"] = {"title": {"text": HEATMAP_UNITS, "side": "right"}}

        module.add_section(
            name=f"Heatmap: {name}",
            anchor=f"tasmanian-heatmap-{slug}",
            description=(
                f"{name} rate at each normalized fragment position for every library, from <code>tasmanian-mismatch</code>. "
                "Rows are libraries and columns are normalized fragment positions."
            ),
            helptext=(
                "Libraries with the same damage pattern have similar rows. Use the plot's cluster switch to group "
                "similar libraries together; positions are never reordered. Positions beyond the length of a "
                "library's longest fragment are empty."
            ),
            plot=plot,
        )
