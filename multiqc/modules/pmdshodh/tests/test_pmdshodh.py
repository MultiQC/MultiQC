import copy

import pytest

from multiqc.modules.pmdshodh.pmdshodh import damage_curves, parse_summary


def parsed(summary):
    d = parse_summary(summary)
    assert d is not None
    return d


def column(values):
    return {b: list(v) for b, v in zip(("A", "C", "G", "T"), values)}


@pytest.fixture
def summary():
    """A version 4 summary with a damage range of 2."""
    zeros = [0, 0]
    return {
        "format_version": 4,
        "tool": {"name": "pmdshodh", "version": "0.3.3"},
        "subcommand": "filter",
        "command": "pmdshodh filter -t 3 --summary s.json s.bam -o k.bam",
        "input": "data/sample1.bam",
        "input_format": "bam",
        "reference": {"path": None, "bases_from": "md", "md_tags": None},
        "mode": {"library": "ds", "treatment": "none", "pmdtools_compat": False},
        "filters": {"threshold": 3.0, "upper_threshold": None},
        "reads": {
            "input": 10,
            "scored": 8,
            "kept": 2,
            "removed_by_threshold": 6,
            "no_informative_sites": 1,
            "skipped": {"unmapped": 2, "too_short": 0},
            "flags": {},
        },
        "score_histogram": {"bin_width": 0.5, "min": -1.0, "max": 0.0, "counts": [6, 2]},
        "length_histogram": {"scored": {"40": 5, "100": 3}, "kept": {"40": 2}},
        "damage": {
            "range": 2,
            "five_prime": {
                "C": column([zeros, [75, 90], zeros, [25, 10]]),
                "G": column([[1, 0], zeros, [99, 100], zeros]),
            },
            "three_prime": {
                "C": column([zeros, [95, 100], zeros, [5, 0]]),
                "G": column([[20, 0], zeros, [80, 0], zeros]),
            },
            "first_position": {"c_to_t": 0.25, "se": 0.04, "n": 100},
        },
    }


def test_counts_and_rates(summary):
    d = parsed(summary)
    assert d["kept"] == 2 and d["removed_by_threshold"] == 6 and d["pct_kept"] == 25.0
    assert d["first_c_to_t"] == 25.0 and d["first_c_to_t_n"] == 100
    assert d["skipped"] == {"unmapped": 2, "too_short": 0}
    assert d["lengths_scored"] == {40: 5, 100: 3} and d["lengths_kept"] == {40: 2}
    # Bin centres, as a percentage of scored reads.
    assert d["score_histogram"] == {-0.75: 75.0, -0.25: 25.0}


def test_damage_ds_uses_g_to_a_at_the_three_prime_end(summary):
    curves = damage_curves(summary)
    assert curves["five_prime"] == {1: 25.0, 2: 10.0}
    # Position 2 at the 3' end has no reference G: no rate rather than 0.
    assert curves["three_prime"] == {1: 20.0}


def test_damage_ss_uses_c_to_t_at_both_ends(summary):
    summary["mode"]["library"] = "ss"
    assert damage_curves(summary)["three_prime"] == {1: 5.0, 2: 0.0}


def test_score_and_damage_runs_have_no_kept_counts(summary):
    summary["subcommand"] = "damage"
    summary["reads"]["kept"] = None
    summary["reads"]["removed_by_threshold"] = None
    summary["length_histogram"]["kept"] = None
    d = parsed(summary)
    assert d["kept"] is None and d["pct_kept"] is None and d["lengths_kept"] == {}


def test_version_1_without_input_format_or_reference(summary):
    old = copy.deepcopy(summary)
    old["format_version"] = 1
    del old["input_format"], old["reference"]
    assert parsed(old)["scored"] == 8


def test_other_tools_and_unknown_versions_are_skipped(summary):
    other = copy.deepcopy(summary)
    other["tool"]["name"] = "pmdtools"
    assert parse_summary(other) is None
    newer = copy.deepcopy(summary)
    newer["format_version"] = 5
    assert parse_summary(newer) is None
