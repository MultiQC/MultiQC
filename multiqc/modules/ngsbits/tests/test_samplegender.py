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
