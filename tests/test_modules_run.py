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


def test_bcftools_stats_multi_sample_hom_het_aggregation(tmp_path, caplog):
    """
    For a multi-sample VCF (one 'set' with several PSC rows), the General Stats hom/het
    values must aggregate across samples rather than keeping only the last PSC row.
    """
    # PSC columns: id sample nRefHom nNonRefHom nHets nTs nTv nIndels depth nSingletons ...
    # nNonRefHom (Hom): 66, 41, 48 ; nHets (Het): 29, 16, 3 ; nRefHom (10/12/15) must be ignored.
    content = """# This file was produced by bcftools stats (1.19+htslib-1.19.1) and can be plotted using plot-vcfstats.
ID\t0\tmulti.vcf.gz
SN\t0\tnumber of records:\t100
SN\t0\tnumber of SNPs:\t90
ST\t0\tA>C\t2
PSC\t0\tS1\t10\t66\t29\t40\t26\t5\t30.0\t3\t0\t0\t0
PSC\t0\tS2\t12\t41\t16\t25\t16\t4\t28.0\t2\t0\t0\t0
PSC\t0\tS3\t15\t48\t3\t20\t10\t2\t25.0\t1\t0\t0\t0
"""

    from multiqc.modules.bcftools.bcftools import MultiqcModule
    from multiqc.types import ColumnKey

    def hom_het(method):
        report.reset()
        config.reset()
        config.bcftools = {"hom_het_aggregation": method} if method else {}
        stats_file = tmp_path / f"multi_{method or 'default'}.stats.txt"
        stats_file.write_text(content)
        report.analysis_files = [stats_file]
        report.search_files(["bcftools"])
        MultiqcModule()
        for section in report.general_stats_data.values():
            for rows in section.values():
                data = rows[0].data
                if ColumnKey("variations_hom") in data:
                    return data[ColumnKey("variations_hom")], data[ColumnKey("variations_het")]
        raise AssertionError("hom/het values not found in General Stats")

    # Default is mean across samples, not the last PSC row (which was 48, 3)
    hom, het = hom_het(None)
    assert hom == 52
    assert het == 16
    assert (hom, het) != (48, 3)

    assert hom_het("mean") == (52, 16)
    assert hom_het("median") == (48, 16)
    assert hom_het("sum") == (155, 48)

    with caplog.at_level("WARNING", logger="multiqc.modules.bcftools.stats"):
        assert hom_het("invalid") == (52, 16)
    assert "Unrecognised bcftools.hom_het_aggregation value 'invalid'" in caplog.text


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
