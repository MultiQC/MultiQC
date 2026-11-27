from pathlib import Path

import pytest

from multiqc import config, report, validation
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.fastdup import MultiqcModule
from multiqc.modules.fastdup.fastdup import parse_report

EXAMPLE = (Path(__file__).parent / "data" / "stats.txt").read_text()
METRICS_ROW = "normal\t169\t49721\t128\t169\t135\t18213\t0\t0.367038\t50018"


@pytest.fixture(autouse=True)
def reset():
    config.reset()
    report.reset()
    validation.reset()
    yield
    report.reset()
    config.reset()
    validation.reset()


def run_module(tmp_path, contents=EXAMPLE, filename="stats.txt"):
    path = tmp_path / filename
    path.write_text(contents)
    config.preserve_module_raw_data = True
    report.analysis_files = [path]
    report.search_files(["fastdup"])
    return MultiqcModule()


def test_reporter_example(tmp_path):
    module = run_module(tmp_path)
    assert module.metrics == {
        "in_test": {
            "LIBRARY": "normal",
            "UNPAIRED_READS_EXAMINED": 169,
            "READ_PAIRS_EXAMINED": 49721,
            "SECONDARY_OR_SUPPLEMENTARY_RDS": 128,
            "UNMAPPED_READS": 169,
            "UNPAIRED_READ_DUPLICATES": 135,
            "READ_PAIR_DUPLICATES": 18213,
            "READ_PAIR_OPTICAL_DUPLICATES": 0,
            "PERCENT_DUPLICATION": 0.367038,
            "ESTIMATED_LIBRARY_SIZE": 50018,
        }
    }
    assert len(module.histograms["in_test"]) == 100
    assert module.histograms["in_test"][1] == 0.999994
    assert module.histograms["in_test"][2] == 1.370063
    assert module.saved_raw_data["multiqc_fastdup"] == module.metrics
    assert module.saved_raw_data["multiqc_fastdup_histogram"] == module.histograms
    assert len(module.sections) == 2
    assert set(report.plot_by_id) == {"fastdup_duplication", "fastdup_yield"}


@pytest.mark.parametrize(
    "command,expected",
    [
        ("fastdup --input /data/alpha.bam --metrics stats.txt", "alpha"),
        ("/usr/local/bin/fastdup --input /data/beta.bam --metrics stats.txt", "beta"),
        ("/home/user/bio tools/fastdup --input /data/install.bam --metrics stats.txt", "install"),
        ("C:\\bin\\FastDup.exe --input C:\\data\\gamma.bam --metrics stats.txt", "gamma"),
        ('fastdup --input "/data/my sample.bam" --metrics stats.txt', "my sample"),
        ("fastdup --input /data/my sample.bam --metrics stats.txt", "stats"),
        ("fastdup --metrics stats.txt --input /data/last.bam", "last"),
        ("fastdup --input=/data/equal.bam --metrics stats.txt", "equal"),
    ],
)
def test_sample_names(tmp_path, command, expected):
    contents = EXAMPLE.replace(EXAMPLE.splitlines()[1], f"# {command}")
    module = run_module(tmp_path, contents)
    assert list(module.metrics) == [expected]


@pytest.mark.parametrize("setting", [True, ["fastdup"]])
def test_use_filename_override(tmp_path, setting):
    config.use_filename_as_sample_name = setting
    assert list(run_module(tmp_path).metrics) == ["stats"]


@pytest.mark.parametrize(
    "setting,expected",
    [
        (False, "in_test"),
        (True, "stats"),
        (["fastdup"], "stats"),
        (["custom_fastdup"], "stats"),
        (["other_module"], "in_test"),
    ],
)
def test_use_filename_override_with_custom_anchor(tmp_path, monkeypatch, setting, expected):
    monkeypatch.setattr(MultiqcModule, "mod_cust_config", {"anchor": "custom_fastdup"})
    config.use_filename_as_sample_name = setting
    module = run_module(tmp_path)
    assert list(module.metrics) == [expected]
    assert list(module.histograms) == [expected]
    assert list(module.saved_raw_data["multiqc_fastdup_custom_fastdup"]) == [expected]
    assert list(module.saved_raw_data["multiqc_fastdup_histogram_custom_fastdup"]) == [expected]


def test_ignore_samples(tmp_path):
    config.sample_names_ignore = ["in_test"]
    with pytest.raises(ModuleNoSamplesFound):
        run_module(tmp_path)


def test_metrics_without_histogram(tmp_path):
    module = run_module(tmp_path, EXAMPLE.split("## HISTOGRAM")[0].rstrip())
    assert module.metrics["in_test"]["READ_PAIR_DUPLICATES"] == 18213
    assert module.histograms["in_test"] == {}
    assert "fastdup_yield" not in report.plot_by_id
    assert module.sections[1].alerts[0].affected_samples == ["in_test"]


def test_zero_reads_and_nonfinite_histogram(tmp_path):
    zero_row = "normal\t0\t0\t0\t0\t0\t0\t0\t0\t0"
    contents = EXAMPLE.split("## HISTOGRAM")[0].replace(METRICS_ROW, zero_row)
    contents += "## HISTOGRAM\tDouble\nBIN CoverageMult\n1\tnan\n2\t-nan\n"
    module = run_module(tmp_path, contents)
    assert module.metrics["in_test"]["PERCENT_DUPLICATION"] == 0
    assert module.histograms["in_test"] == {1: None, 2: None}
    assert module.saved_raw_data["multiqc_fastdup_histogram"]["in_test"] == {1: None, 2: None}
    assert not report.plot_by_id
    assert module.sections[0].alerts[0].affected_samples == ["in_test"]
    assert module.sections[1].alerts


def test_unestimable_library_does_not_plot_zero_yield(tmp_path):
    contents = EXAMPLE.replace(METRICS_ROW, METRICS_ROW.rsplit("\t", 1)[0] + "\t0")
    module = run_module(tmp_path, contents)
    assert len(module.histograms["in_test"]) == 100
    assert "fastdup_yield" not in report.plot_by_id
    assert "fastdup_duplication" in report.plot_by_id


def test_partial_nonfinite_histogram(tmp_path):
    contents = EXAMPLE.replace("1\t0.999994", "1\tnan")
    module = run_module(tmp_path, contents)
    assert module.histograms["in_test"][1] is None
    assert module.histograms["in_test"][2] == 1.370063
    assert "fastdup_yield" in report.plot_by_id
    assert module.sections[1].alerts[-1].affected_samples == ["in_test"]


def test_preserve_numeric_read_group_and_large_integer():
    large = 2**53 + 1
    row = f"001\t{large}\t0\t0\t0\t1\t0\t0\t0\t0"
    metrics, _ = parse_report(EXAMPLE.replace(METRICS_ROW, row))
    assert metrics["LIBRARY"] == "001"
    assert metrics["UNPAIRED_READS_EXAMINED"] == large


@pytest.mark.parametrize(
    "contents",
    [
        EXAMPLE.replace(METRICS_ROW, "normal\t169"),
        EXAMPLE.replace("READ_PAIRS_EXAMINED", "MISSING_PAIRS"),
        EXAMPLE.replace("\t0.367038\t", "\tnan\t"),
        EXAMPLE.replace("\t0.367038\t", "\t1.1\t"),
        EXAMPLE.replace(METRICS_ROW, METRICS_ROW.replace("\t169\t", "\t-1\t", 1)),
        EXAMPLE.replace("\t18213\t0\t", "\t18213\t18214\t"),
        EXAMPLE.replace(METRICS_ROW, METRICS_ROW + "\n" + METRICS_ROW),
        "# fastdup --input input.bam\n## METRICS\nLIBRARY\tPERCENT_DUPLICATION\n",
        EXAMPLE.split("## HISTOGRAM")[0] + "## HISTOGRAM\tDouble\n",
        EXAMPLE.replace("BIN CoverageMult", "BIN WrongColumn"),
        EXAMPLE.replace("1\t0.999994", "1\t0.999994\t12"),
        EXAMPLE.replace("1\t0.999994", "1\t0.999994\n1\t0.999994"),
    ],
)
def test_malformed_report_fails_with_file_context(tmp_path, contents):
    with pytest.raises(ValueError, match="Invalid FastDup report .*stats.txt:"):
        run_module(tmp_path, contents)


@pytest.mark.parametrize("header", ["# picard.sam.MarkDuplicates INPUT=x.bam", "# notfastdup --input x.bam"])
def test_unrelated_reports_not_detected(tmp_path, header):
    contents = EXAMPLE.replace(EXAMPLE.splitlines()[1], header)
    with pytest.raises(ModuleNoSamplesFound):
        run_module(tmp_path, contents)


def test_duplicate_sample_overwrites_metrics_and_histogram_together(tmp_path):
    first = tmp_path / "first.txt"
    second = tmp_path / "second.txt"
    first.write_text(EXAMPLE)
    second.write_text(EXAMPLE.split("## HISTOGRAM")[0].replace("\t0.367038\t", "\t0.4\t"))
    config.preserve_module_raw_data = True
    report.analysis_files = [first, second]
    report.search_files(["fastdup"])
    module = MultiqcModule()
    assert len(module.metrics) == 1
    assert module.metrics["in_test"]["PERCENT_DUPLICATION"] == 0.4
    assert module.histograms["in_test"] == {}


def test_optical_duplicates_are_not_double_counted(tmp_path):
    module = run_module(tmp_path, EXAMPLE.replace("\t18213\t0\t", "\t18213\t101\t"))
    plot = report.plot_by_id["fastdup_duplication"].model_dump()
    counts = {category["name"]: category["data"][0] for category in plot["datasets"][0]["cats"]}
    assert counts == {
        "Unique paired reads": 63016,
        "Unique unpaired reads": 34,
        "Optical duplicate paired reads": 202,
        "Other duplicate paired reads": 36224,
        "Duplicate unpaired reads": 135,
        "Unmapped reads": 169,
    }
    assert sum(counts.values()) == 99780
    assert module.metrics["in_test"]["READ_PAIR_DUPLICATES"] == 18213
    assert module.metrics["in_test"]["READ_PAIR_OPTICAL_DUPLICATES"] == 101
    line = report.plot_by_id["fastdup_yield"].model_dump()["datasets"][0]["lines"][0]
    assert line["pairs"] == list(module.histograms["in_test"].items())


def test_custom_anchor_scopes_plot_ids_and_exports(tmp_path, monkeypatch):
    monkeypatch.setattr(MultiqcModule, "mod_cust_config", {"anchor": "custom_fastdup"})
    module = run_module(tmp_path)
    assert set(report.plot_by_id) == {"custom_fastdup_duplication", "custom_fastdup_yield"}
    assert set(module.saved_raw_data) == {
        "multiqc_fastdup_custom_fastdup",
        "multiqc_fastdup_histogram_custom_fastdup",
    }
