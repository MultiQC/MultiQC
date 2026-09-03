import tempfile
from pathlib import Path
from typing import Callable, List, Union

import pytest

from multiqc import BaseMultiqcModule, config, parse_logs, report, reset
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.core.update_config import ClConfig, update_config
from multiqc.types import SectionKey

modules = [(k, entry_point) for k, entry_point in config.avail_modules.items() if k != "custom_content"]


@pytest.fixture(scope="module")
def multiqc_reset():
    reset()


@pytest.fixture(autouse=True)
def reset_config():
    """Reset config state after each test."""
    original_strict = config.strict
    original_sample_names_ignore = config.sample_names_ignore[:]
    yield
    config.strict = original_strict
    config.sample_names_ignore[:] = original_sample_names_ignore


@pytest.mark.parametrize("module_id,entry_point", modules)
def test_all_modules(module_id, entry_point, data_dir):
    """
    Verify that all modules do at least something.

    Modules can add to this test by providing more specific actions for specific files:
    e.g. if a file should be skipped, or cause a runtime error, or the raw data should
    look in a specific way.
    """

    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()

    config.strict = True

    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    _module = module_cls()
    for m in _module if isinstance(_module, List) else [_module]:
        assert len(report.general_stats_data) > 0 or len(m.sections) > 0


def test_bcftools_stats_zero_depth_samples(tmp_path):
    """All-zero average depth values should still render a sequencing-depth plot."""
    stats_file = tmp_path / "rotavirus.stats.txt"
    stats_file.write_text(
        """# This file was produced by bcftools stats (1.19+htslib-1.19.1) and can be plotted using plot-vcfstats.
ID\t0\trotavirus.vcf.gz
SN\t0\tnumber of records:\t2
SN\t0\tnumber of SNPs:\t2
ST\t0\tA>C\t2
PSC\t0\tS1\t0\t0\t1\t1\t1\t0\t0.0\t0\t0\t0\t0
PSC\t0\tS2\t0\t0\t1\t1\t1\t0\t0.0\t0\t0\t0\t0
"""
    )

    report.analysis_files = [stats_file]
    report.search_files(["bcftools"])

    from multiqc.modules.bcftools.bcftools import MultiqcModule

    module = MultiqcModule()

    assert "bcftools-stats_sequencing_depth" in {section.id for section in module.sections}
    assert "bcftools-stats-sequencing-depth" in report.plot_by_id


@pytest.mark.parametrize("module_id,entry_point", modules)
def test_ignore_samples(module_id, entry_point, data_dir):
    """
    Verify all modules call self.ignore_samples() correctly.
    """
    mod_dir = data_dir / "modules" / module_id
    assert mod_dir.exists() and mod_dir.is_dir()

    config.sample_names_ignore = ["*"]

    report.analysis_files = [mod_dir]
    report.search_files([module_id])

    module_cls: Callable[[], Union[BaseMultiqcModule, List[BaseMultiqcModule]]] = entry_point.load()
    with pytest.raises(ModuleNoSamplesFound):
        _module = module_cls()


@pytest.mark.parametrize(
    ["config_options", "expected_to_write"],
    [
        (dict(), True),
        ({"make_data_dir": False}, False),
        ({"filename": "stdout"}, False),
    ],
)
def test_write_data_file(monkeypatch, tmp_path, config_options, expected_to_write):
    """
    Test module.write_data_file() write something
    """
    (tmp_path / "multiqc_tmp").mkdir()
    monkeypatch.setattr(tempfile, "mkdtemp", lambda: tmp_path / "multiqc_tmp")
    config.update(config_options)
    module = BaseMultiqcModule()
    module.write_data_file({"Sample": {"key": "value"}}, "multiqc_mymodule")

    expected_path = tmp_path / "multiqc_tmp" / "multiqc_data" / "multiqc_mymodule.txt"
    if expected_to_write:
        assert expected_path.exists()
        assert expected_path.open().read().strip() == """Sample\tkey\nSample\tvalue""".strip()
    else:
        assert not expected_path.exists()


@pytest.mark.parametrize(
    "use_filename_as_sample_name,fn_clean_sample_names,prepend_dirs,expected_sample_name",
    [
        (None, None, None, "SAMPLE_FROM_CONTENTS"),
        (True, None, None, "SAMPLE_FROM_FILENAME"),
        (None, False, None, "SAMPLE_FROM_CONTENTS.fastq.gz"),
        (True, False, None, "SAMPLE_FROM_FILENAME.stderr"),
        (None, None, True, "subdir | SAMPLE_FROM_CONTENTS"),
        (True, None, True, "subdir | SAMPLE_FROM_FILENAME"),
        (["trimmomatic"], None, None, "SAMPLE_FROM_FILENAME"),
        (["other_module"], None, None, "SAMPLE_FROM_CONTENTS"),  # Should not affect trimmomatic
    ],
)
def test_use_filename_as_sample_name(
    multiqc_reset, tmp_path, use_filename_as_sample_name, fn_clean_sample_names, prepend_dirs, expected_sample_name
):
    """
    Verify that `--fn_as_s_name`, `--fullnames`, and `--dirs` works
    """
    report.reset()

    MODULE_NAME = "trimmomatic"
    (tmp_path / "subdir").mkdir()
    input_file = tmp_path / "subdir" / "SAMPLE_FROM_FILENAME.stderr"
    input_file.write_text("""\
TrimmomaticSE: Started with arguments: SAMPLE_FROM_CONTENTS.fastq.gz
Input Reads: 39733090 Surviving: 32590558 (82.02%) Dropped: 7142532 (17.98%)
TrimmomaticSE: Completed successfully""")

    update_config(
        cfg=ClConfig(
            run_modules=[MODULE_NAME],
            use_filename_as_sample_name=use_filename_as_sample_name,
            fn_clean_sample_names=fn_clean_sample_names,
            prepend_dirs=prepend_dirs,
            dirs_depth=1 if prepend_dirs else None,
            preserve_module_raw_data=True,
        )
    )

    report.analysis_files = [input_file]
    report.search_files([MODULE_NAME])

    from multiqc.modules.trimmomatic.trimmomatic import MultiqcModule

    m = MultiqcModule()

    assert m.saved_raw_data is not None
    assert expected_sample_name in m.saved_raw_data[f"multiqc_{MODULE_NAME}"]


@pytest.mark.parametrize(
    # Custom anchors are used to suffix the write_data_file fn (= saved_raw_data key)
    # If custom anchor is not provided for a repeated module, it is added for the second occurent
    # by sanitising the module id (e.g. adapterremoval -> adapterremoval-1
    "anchors,expected_raw_data_keys",
    [
        (
            ["my_anchor_se", "my_anchor_pe"],
            ["multiqc_adapter_removal_my_anchor_se", "multiqc_adapter_removal_my_anchor_pe"],
        ),
        ([None, None], ["multiqc_adapter_removal", "multiqc_adapter_removal_adapterremoval-1"]),
    ],
)
def test_path_filters(multiqc_reset, tmp_path, data_dir, anchors, expected_raw_data_keys):
    search_path = data_dir / "modules" / "adapterremoval"
    assert search_path.exists() and search_path.is_dir()

    expected_pe_files = {
        "paired_end_collapsed/pec1.settings",
        "paired_end_collapsed/pec2.settings",
        "paired_end_noncollapsed/penc1.settings",
        "paired_end_noncollapsed/penc2.settings",
    }
    expected_se_files = {
        "single_end/se.settings",
    }

    assert all((search_path / fn).exists() for fn in expected_pe_files)
    assert all((search_path / fn).exists() for fn in expected_pe_files)

    parse_logs(
        search_path,
        module_order=[
            {
                "adapterremoval": {
                    "name": "adapterremoval (single end)",
                    "anchor": anchors[0],
                    "path_filters": ["*/se.*"],
                },
            },
            {
                "adapterremoval": {
                    "name": "adapterremoval (paired end)",
                    "anchor": anchors[1],
                    "path_filters": ["*/pec?.*", "*/penc?.*"],
                },
            },
        ],
        preserve_module_raw_data=True,
        strict=True,
    )

    assert len(report.modules) == 2
    assert len(report.general_stats_data) == 2
    assert report.modules[0].name == "adapterremoval (single end)"
    assert report.modules[1].name == "adapterremoval (paired end)"
    assert report.modules[0].saved_raw_data is not None
    assert report.modules[1].saved_raw_data is not None
    assert report.modules[0].saved_raw_data[expected_raw_data_keys[0]].keys() == {
        Path(fn).name for fn in expected_se_files
    }
    assert report.modules[1].saved_raw_data[expected_raw_data_keys[1]].keys() == {
        Path(fn).name for fn in expected_pe_files
    }
    assert set(report.general_stats_data[SectionKey(anchors[0] or "adapterremoval")].keys()) == {
        Path(fn).name for fn in expected_se_files
    }
    assert set(report.general_stats_data[SectionKey(anchors[1] or "adapterremoval-1")].keys()) == {
        Path(fn).name for fn in expected_pe_files
    }


def _run_freyja(tmp_path):
    from multiqc.modules.freyja.freyja import MultiqcModule

    report.reset()
    report.analysis_files = [tmp_path]
    report.search_files(["freyja"])
    return MultiqcModule()


def _freyja_stats_cols(sample):
    """Column keys recorded across all freyja general-stats sections for a given sample."""
    cols: set = set()
    for section_key, section in report.general_stats_data.items():
        if not str(section_key).startswith("freyja"):
            continue
        for row in section.get(sample, []):
            cols.update(row.data.keys())
    return cols


def test_freyja_legacy_no_pathogen(tmp_path):
    """Legacy Freyja reports without a 'pathogen' line default to SARS-CoV-2."""
    (tmp_path / "sample.demix_out.tsv").write_text(
        "\tvariant_files/sample.demix_out.tsv\n"
        "summarized\t[('BA.5* [Omicron (BA.5.X)]', 0.9), ('Omicron', 0.1)]\n"
        "lineages\tBF.5.1 BF.7.7\n"
        "abundances\t0.5 0.5\n"
        "resid\t0.06\n"
        "coverage\t1.83\n"
    )

    module = _run_freyja(tmp_path)

    # The pathogen (defaulted to SARS-CoV-2) is included in anchors and column keys
    assert {s.anchor for s in module.sections} == {"freyja-summary-sars-cov-2"}
    assert _freyja_stats_cols("sample.demix_out") == {
        "Top_lineage_freyja_sars_cov_2",
        "Top_lineage_freyja_sars_cov_2_percentage",
    }


def test_freyja_multiple_pathogens(tmp_path):
    """When reports declare different pathogens, each gets its own section and columns."""
    (tmp_path / "covid.demix_out.tsv").write_text(
        "\tvariant_files/covid.demix_out.tsv\n"
        "summarized\t[('BA.5* [Omicron (BA.5.X)]', 0.9), ('Omicron', 0.1)]\n"
        "lineages\tBF.5.1 BF.7.7\n"
        "abundances\t0.5 0.5\n"
        "resid\t0.06\n"
        "coverage\t1.83\n"
        "pathogen\tSARS-CoV-2\n"
    )
    # A legacy report with no pathogen line should be grouped with SARS-CoV-2
    (tmp_path / "covid_legacy.demix_out.tsv").write_text(
        "\tvariant_files/covid_legacy.demix_out.tsv\n"
        "summarized\t[('BA.5* [Omicron (BA.5.X)]', 0.8), ('Omicron', 0.2)]\n"
        "lineages\tBF.5.1\n"
        "abundances\t1.0\n"
        "resid\t0.1\n"
        "coverage\t1.83\n"
    )
    (tmp_path / "rsv.demix_out.tsv").write_text(
        "\tvariant_files/rsv.demix_out.tsv\n"
        "summarized\t[('A.D.1', 0.7), ('A.D.5.2', 0.3)]\n"
        "lineages\tA.D.1 A.D.5.2\n"
        "abundances\t0.7 0.3\n"
        "resid\t0.05\n"
        "coverage\t1.9\n"
        "pathogen\tRSV-A\n"
    )

    module = _run_freyja(tmp_path)

    assert {s.anchor for s in module.sections} == {"freyja-summary-sars-cov-2", "freyja-summary-rsv-a"}

    # The legacy sample (no pathogen line) is grouped under SARS-CoV-2 alongside the explicit one
    assert _freyja_stats_cols("covid.demix_out") == {
        "Top_lineage_freyja_sars_cov_2",
        "Top_lineage_freyja_sars_cov_2_percentage",
    }
    assert _freyja_stats_cols("covid_legacy.demix_out") == {
        "Top_lineage_freyja_sars_cov_2",
        "Top_lineage_freyja_sars_cov_2_percentage",
    }
    assert _freyja_stats_cols("rsv.demix_out") == {
        "Top_lineage_freyja_rsv_a",
        "Top_lineage_freyja_rsv_a_percentage",
    }
