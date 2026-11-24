from pathlib import Path

import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.fgumi import MultiqcModule

FAMILY_SIZES_TSV = (
    "family_size\tcount\tfraction\tfraction_gt_or_eq_family_size\n1\t60\t0.6\t1\n2\t30\t0.3\t0.4\n5\t10\t0.1\t0.1\n"
)


def test_family_sizes_parsed_and_plotted(run_fgumi):
    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1"}
    assert module.saved_raw_data["multiqc_fgumi_family_sizes"]["S1"] == {1: 60, 2: 30, 5: 10}
    assert "fgumi_family_sizes" in report.plot_by_id


def test_user_named_histogram_is_detected_by_header(run_fgumi):
    module = run_fgumi({"my_hist.tsv": FAMILY_SIZES_TSV})
    assert module.samples_parsed_by_tool["family_sizes"] == {"my_hist"}


@pytest.mark.parametrize(
    "fn",
    [
        "S1.family_sizes.txt",
        "S1.duplex_family_sizes.txt",  # a longer suffix must not leave "S1.duplex_family_sizes"
        "S1.grouping_metrics.txt",
        "S1.umi_counts.txt",
        # `fgumi runall --all-metrics S1` puts the stage before the suffix
        "S1.group.family_sizes.txt",
        "S1.duplex.umi_counts.txt",
        "S1.simplex.simplex_yield_metrics.txt",
        "S1.codec.family_sizes.txt",
        "S1.correct.metrics.txt",
        "S1.filter.stats.txt",
    ],
)
def test_sample_name_strips_fgumi_suffixes(run_fgumi, fn):
    from multiqc.modules.fgumi.util import sample_name

    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    f = {"fn": fn, "root": "", "s_name": fn, "sp_key": "fgumi/family_sizes"}  # the keys find_log_files provides
    assert sample_name(module, f) == "S1"


def test_fullnames_keeps_the_file_name(run_fgumi):
    from multiqc import config

    config.fn_clean_sample_names = False
    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES_TSV})
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1.family_sizes.txt"}


def test_runall_outputs_share_one_sample(run_fgumi):
    grouping = "\t".join(
        ["accepted_sam_records", "discarded_non_pf", "discarded_poor_alignment", "discarded_ns_in_umi"]
        + ["discarded_umis_to_short"]
    )
    module = run_fgumi(
        {
            "S1.group.family_sizes.txt": FAMILY_SIZES_TSV,
            "S1.group.grouping_metrics.txt": grouping + "\n10\t0\t0\t0\t0\n",
        }
    )
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1"}
    assert module.samples_parsed_by_tool["grouping"] == {"S1"}
    # Each file is listed as a source for the sample; neither replaces the other.
    sources = {s for s, by_sample in report.data_sources["fgumi"].items() if "S1" in by_sample}
    assert {"family_sizes", "grouping_metrics"} <= sources


FGBIO_ERROR_RATE = (
    "read_number\tposition\tbases_total\terrors\terror_rate\ta_to_c_error_rate\ta_to_g_error_rate"
    "\ta_to_t_error_rate\tc_to_a_error_rate\tc_to_g_error_rate\tc_to_t_error_rate\n"
    "1\t1\t1000\t0\t0.0\t0.001\t0.002\t0.003\t0.004\t0.005\t0.006\n"
)
POSITION_GROUP_SIZES = "position_group_size\tcount\tfraction\tfraction_gt_or_eq_position_group_size\n1\t10\t1\t1\n"


HIST = FAMILY_SIZES_TSV
FGUMI_ONLY = POSITION_GROUP_SIZES  # a file only fgumi writes
FGBIO_ONLY = FGBIO_ERROR_RATE  # a file only fgbio writes


@pytest.mark.parametrize(
    "layout, family_sizes_module, modules, expected",
    [
        pytest.param({"A.family_sizes.txt": HIST}, None, ["fgbio", "fgumi"], {"A": "fgbio"}, id="no-evidence-fgbio"),
        pytest.param(
            {"A.family_sizes.txt": HIST, "A.position_group_sizes.txt": FGUMI_ONLY},
            None,
            ["fgbio", "fgumi"],
            {"A": "fgumi"},
            id="same-dir-fgumi-evidence",
        ),
        pytest.param(
            {"A.family_sizes.txt": HIST, "A.error_rate.txt": FGBIO_ONLY},
            None,
            ["fgbio", "fgumi"],
            {"A": "fgbio"},
            id="same-dir-fgbio-evidence",
        ),
        pytest.param(
            {"A.family_sizes.txt": HIST, "A.position_group_sizes.txt": FGUMI_ONLY, "A.error_rate.txt": FGBIO_ONLY},
            None,
            ["fgbio", "fgumi"],
            {"A": "fgbio"},
            id="mixed-dir-fgbio",
        ),
        pytest.param(
            {"S1/group/A.family_sizes.txt": HIST, "S1/dedup/A.position_group_sizes.txt": FGUMI_ONLY},
            None,
            ["fgbio", "fgumi"],
            {"A": "fgumi"},
            id="evidence-in-parent-dir",
        ),
        pytest.param(
            {
                "fgbio/B.family_sizes.txt": HIST,
                "fgbio/B.error_rate.txt": FGBIO_ONLY,
                "fgumi/A.family_sizes.txt": HIST,
                "fgumi/A.position_group_sizes.txt": FGUMI_ONLY,
            },
            None,
            ["fgbio", "fgumi"],
            {"A": "fgumi", "B": "fgbio"},
            id="sibling-dirs-each-keep-their-own",
        ),
        pytest.param(
            {
                "A.error_rate.txt": FGBIO_ONLY,
                "fgumi/A.family_sizes.txt": HIST,
                "fgumi/A.position_group_sizes.txt": FGUMI_ONLY,
            },
            None,
            ["fgbio", "fgumi"],
            {"A": "fgumi"},
            id="nearest-dir-wins-over-mixed-parent",
        ),
        pytest.param({"A.family_sizes.txt": HIST}, "fgumi", ["fgbio", "fgumi"], {"A": "fgumi"}, id="config-fgumi"),
        pytest.param(
            {"A.family_sizes.txt": HIST, "A.position_group_sizes.txt": FGUMI_ONLY},
            "fgbio",
            ["fgbio", "fgumi"],
            {"A": "fgbio"},
            id="config-fgbio",
        ),
        pytest.param({"A.family_sizes.txt": HIST}, None, ["fgumi"], {"A": "fgumi"}, id="only-fgumi-running"),
        pytest.param({"A.family_sizes.txt": HIST}, "fgumi", ["fgbio"], {"A": "fgbio"}, id="only-fgbio-running"),
    ],
)
def test_one_module_reports_each_family_size_histogram(tmp_path, layout, family_sizes_module, modules, expected):
    # fgumi and fgbio GroupReadsByUmi histograms are identical; exactly one module may report each file.
    from multiqc import config
    from multiqc.modules.fgbio import MultiqcModule as FgbioModule

    config.reset()
    if family_sizes_module is not None:
        config.fgumi_config = {"family_sizes_module": family_sizes_module}
    paths = []
    for name, content in layout.items():
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        paths.append(path)
    report.reset()
    report.analysis_files = paths
    report.search_files(modules)
    try:
        for module_class in [FgbioModule, MultiqcModule]:
            if module_class.__module__.split(".")[2] in modules:
                try:
                    module_class()
                except ModuleNoSamplesFound:
                    pass
        # The module that recorded each histogram as a data source is the one that reported it.
        owners = {
            Path(path).name.split(".")[0]: module_name
            for module_name, sections in report.data_sources.items()
            for by_sample in sections.values()
            for path in by_sample.values()
            if str(path).endswith(".family_sizes.txt")
        }
        assert owners == expected
    finally:
        config.reset()


def test_malformed_file_is_skipped(run_fgumi, caplog):
    truncated = FAMILY_SIZES_TSV + "7\t3\n"
    with pytest.raises(ModuleNoSamplesFound):
        run_fgumi({"S1.family_sizes.txt": truncated})
    assert any("S1.family_sizes.txt" in record.message for record in caplog.records)


def test_no_fgumi_files_raises(run_fgumi):
    with pytest.raises(ModuleNoSamplesFound):
        run_fgumi({"irrelevant.txt": "hello world\n"})
