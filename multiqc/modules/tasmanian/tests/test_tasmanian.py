import pytest

from multiqc.modules.tasmanian.inconsistencies import parse_inconsistencies_table
from multiqc.modules.tasmanian.mismatch import (
    ALL_MISMATCHES,
    _profile_rates,
    class_rate,
    overall_rates,
    parse_mismatch_table,
)
from multiqc.modules.tasmanian.variants import parse_variants_table

RAW = """base_change\tread_num\treference_order\tfragment_position\tcount
C>C\t1\t1\t1\t90
C>T\t1\t1\t1\t10
G>G\t1\t1\t1\t95
G>T\t1\t1\t1\t5
C>C\t1\t2\t1\t10
C>T\t1\t2\t1\t10
C>C\t2\t2\t1\t100
"""


def test_parse_raw_sums_reference_order():
    data = parse_mismatch_table(RAW)
    assert data["profile"][1][1] == {"C>C": 100, "C>T": 20, "G>G": 95, "G>T": 5}
    assert data["profile"][2][1] == {"C>C": 100}


def test_parse_header_only():
    data = parse_mismatch_table("base_change\tread_num\treference_order\tfragment_position\tcount\n")
    assert data["profile"] == {}


@pytest.mark.parametrize(
    "table",
    [
        "base_change\tread_num\treference_order\tread_position\tcount\nC>T\t1\t1\t1\t5\n",
        "base_change\tread_num\treference_order\tfragment_position\tnormalized_frequency\nC>T\t1\t1\t1\t0.5\n",
        "base_change\tread_num\treference_order\tfragment_position\tcount\nC>N\t1\t1\t1\t5\n",
        "base_change\tread_num\treference_order\tfragment_position\tcount\nC>T\t1\t1\t1\n",
        "base_change\tread_num\treference_order\tfragment_position\tcount\nC>T\t1\t1\t1\tmany\n",
    ],
)
def test_parse_mismatch_invalid(table):
    with pytest.raises(ValueError):
        parse_mismatch_table(table)


def test_overall_rates():
    rates = overall_rates(parse_mismatch_table(RAW))
    assert rates is not None
    # 25 of 320 bases mismatch
    assert rates[ALL_MISMATCHES] == pytest.approx(100 * 25 / 320)
    assert rates["total_bases"] == 320
    # C>T: 20 of 220 bases with reference C
    assert rates["C>T"] == pytest.approx(100 * 20 / 220)
    assert rates["G>T"] == pytest.approx(5.0)
    # No reference A in the table
    assert "A>C" not in rates


def test_class_rate_pools_reference_bases():
    counts = {"C>C": 90, "C>T": 10, "G>G": 45, "G>A": 5, "G>T": 5, "A>A": 1000}
    # Pooled over reference C (100 bases) and G (55 bases), not an average of 10% and 9.1%
    assert class_rate(counts, ["C>T", "G>A"]) == pytest.approx(100 * 15 / 155)
    assert class_rate(counts, ["G>T", "C>A"]) == pytest.approx(100 * 5 / 155)
    assert class_rate(counts, ["C>T"]) == pytest.approx(10.0)
    # No bases with reference T
    assert class_rate(counts, ["T>C"]) is None


def test_overall_rates_include_groups():
    rates = overall_rates(parse_mismatch_table(RAW))
    assert rates is not None
    # C>T 20 and G>A 0 over reference C (220) and G (100)
    assert rates["Deamination (C>T + G>A)"] == pytest.approx(100 * 20 / 320)
    # G>T 5 and C>A 0 over the same bases
    assert rates["Oxidation (G>T + C>A)"] == pytest.approx(100 * 5 / 320)


def test_profile_rates_groups():
    rates = _profile_rates(parse_mismatch_table(RAW))
    assert rates["Deamination (C>T + G>A)"] == {1: pytest.approx(100 * 20 / 320)}


def test_profile_rates_merges_reads():
    data = parse_mismatch_table(
        "base_change\tread_num\treference_order\tfragment_position\tcount\n"
        "C>C\t1\t1\t5\t90\nC>T\t1\t1\t5\t10\nC>C\t2\t2\t5\t70\nC>T\t2\t2\t5\t30\n"
    )
    assert _profile_rates(data)["C>T"] == {5: pytest.approx(20.0)}


def test_parse_variants():
    data = parse_variants_table(
        "chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth\n"
        "chr1\t10\tC\tT\t12\t38\nchr1\t20\tC\tT\t9\t27\nchr2\t30\tG\tA\t15\t41\n"
    )
    assert data["variant_sites"] == 3
    assert data["mismatch_count"] == 36
    assert data["sites_by_class"]["C>T"] == 2
    assert data["sites_by_class"]["G>A"] == 1


def test_parse_variants_header_only():
    data = parse_variants_table("chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth\n")
    assert data["variant_sites"] == 0
    assert data["mismatch_count"] == 0


def test_parse_variants_invalid_substitution():
    with pytest.raises(ValueError):
        parse_variants_table("chromosome\tposition\treference_base\tmismatch_base\tcount\tdepth\nchr1\t1\tC\tC\t1\t1\n")


def test_parse_inconsistencies():
    data = parse_inconsistencies_table(
        "read1_position\tread2_position\tdiscordance_type\tcount\n"
        "18\t83\tR1:A_R2:G\t5\n19\t82\tR1:C_R2:T\t3\n20\t81\tR1:A_R2:G\t2\n"
    )
    assert data["inconsistent_bases"] == 10
    assert data["by_type"] == {"R1:A_R2:G": 7, "R1:C_R2:T": 3}


def test_parse_inconsistencies_invalid_type():
    with pytest.raises(ValueError):
        parse_inconsistencies_table("read1_position\tread2_position\tdiscordance_type\tcount\n1\t2\t<b>\t1\n")
