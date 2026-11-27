import pytest

from multiqc.modules.ngsbits.samplegender import parse_file


def test_parse_samplegender_report():
    report = """\
#file\tgender\treads_chry\treads_chrx\tratio_chry_chrx
SAMPLE1_sorted_md.bam\tmale\t3872086\t14244970\t0.2718
"""

    assert parse_file(report) == {
        "gender": "male",
        "reads_chry": 3872086.0,
        "reads_chrx": 14244970.0,
        "ratio_chry_chrx": 0.2718,
    }


def test_parse_samplegender_unknown_gender_with_tabbed_description():
    report = """\
#sample\tgender\treads_chry\treads_chrx\tratio_chry_chrx
XXX\tunknown\t(ratio\tin\tgray\tarea)\t12246\t153042\t0.0800
"""

    assert parse_file(report) == {
        "gender": "unknown (ratio in gray area)",
        "reads_chry": 12246.0,
        "reads_chrx": 153042.0,
        "ratio_chry_chrx": pytest.approx(0.08),
    }


def test_parse_samplegender_unknown_gender_with_spaces():
    report = """\
#file\tgender\treads_chry\treads_chrx\tratio_chry_chrx
S.bam\tunknown (ratio in gray area)\t12246\t153042\t0.0800
"""

    assert parse_file(report) == {
        "gender": "unknown (ratio in gray area)",
        "reads_chry": 12246.0,
        "reads_chrx": 153042.0,
        "ratio_chry_chrx": pytest.approx(0.08),
    }


def test_parse_samplegender_hetx_method():
    report = """\
#file\tgender\tsnps_usable\thom_count\thet_count\thet_fraction
S.bam\tmale\t123 of 456\t100\t23\t0.1870
"""

    assert parse_file(report) == {
        "gender": "male",
        "snps_usable": "123 of 456",
        "hom_count": 100.0,
        "het_count": 23.0,
        "het_fraction": pytest.approx(0.187),
    }


def test_parse_samplegender_truncated_row():
    report = """\
#file\tgender\treads_chry\treads_chrx\tratio_chry_chrx
S\tmale\t1\t2
"""

    with pytest.raises(ValueError):
        parse_file(report)
