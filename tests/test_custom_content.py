import tempfile
from pathlib import Path

import pytest

import multiqc
from multiqc import report, config
from multiqc.core.update_config import update_config, ClConfig
from multiqc.modules.custom_content import custom_module_classes
from multiqc.utils import testing
from multiqc.validation import ConfigValidationError
from multiqc.core.file_search import file_search


def test_custom_content(tmp_path):
    file = tmp_path / "mysample_mqc.txt"
    id = "dupradar"
    file.write_text(
        f"""\
#id: {id}
#plot_type: 'linegraph'
#section_name: 'DupRadar'
#section_href: 'bioconductor.org/packages/release/bioc/html/dupRadar.html'
#description: "first line.
#    Second line "
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xlog: True
#    xlab: 'expression (reads/kbp)'
#    ymax: 100
#    ymin: 0
#    tt_label: '<b>{{point.x:.1f}} reads/kbp</b>: {{point.y:,.2f}}% duplicates'
#    x_line:
#        - color: 'green'
#          dash: 'longdash'
#          label: '0.5 RPKM'
#          value: 0.5
#          width: 1
#    x_bands:
#        - from: 0.1
#          to: 0.3
#          color: 'black'
#    y_bands:
#        - from: 0.1
#          to: 0.3
0.561167227833894 0.0146313784854042
3.63901018922853 0.0639394516274346
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert f"{id}-section-plot" in report.plot_by_id
    assert report.plot_by_id[f"{id}-section-plot"].id == f"{id}-section-plot"
    assert report.plot_by_id[f"{id}-section-plot"].plot_type == "xy_line"


def test_deprecated_fields(tmp_path, capsys):
    file = tmp_path / "mysample_mqc.txt"
    id = "dupradar"
    file.write_text(
        f"""\
#id: {id}
#plot_type: 'linegraph'
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xLog: True
#    xPlotLines:
#        - colour: 'green'
#          dash: 'LongDash'
#          label:
#                style: {{color: 'green'}}
#                text: '0.5 RPKM'
#                verticalAlign: 'bottom'
#                y: -65
#          value: 0.5
#          width: 1
0.561167227833894 0.0146313784854042
3.63901018922853 0.0639394516274346
"""
    )

    report.analysis_files = [file]
    update_config(cfg=ClConfig(verbose=True))  # force re-creating logger
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert f"{id}-section-plot" in report.plot_by_id
    assert report.plot_by_id[f"{id}-section-plot"].id == f"{id}-section-plot"
    assert report.plot_by_id[f"{id}-section-plot"].plot_type == "xy_line"

    err = str(capsys.readouterr().err)
    assert "Line plot's x_lines or y_lines 'label' field is expected to be a string" in err
    assert "'LongDash' is a deprecated dash style, use 'longdash'" in err
    assert "Deprecated field 'colour'. Use 'color' instead" in err
    assert "Deprecated field 'xLog'. Use 'xlog' instead" in err
    assert "Deprecated field 'xPlotLines'. Use 'x_lines' instead" in err


@pytest.mark.parametrize("strict", [True, False])
def test_wrong_fields(tmp_path, capsys, strict, monkeypatch):
    """
    Values of wrong types. Should fail in strict mode, but still produce output in non-strict mode.
    """
    file = tmp_path / "mysample_mqc.txt"
    id = "dupradar"
    file.write_text(
        f"""\
#id: {id}
#plot_type: 'linegraph'
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xLog: True
#    xlab: True
#    y__lab: '% duplicate reads'
#    ymax: 100
#    ymin: "0"
0.561167227833894 0.0146313784854042
"""
    )
    (tmp_path / "tmp").mkdir()
    monkeypatch.setattr(tempfile, "mkdtemp", lambda: tmp_path / "tmp")

    config.strict = strict
    report.analysis_files = [file]
    report.search_files(["custom_content"])

    if strict:
        with pytest.raises(ConfigValidationError):
            custom_module_classes()
    else:
        custom_module_classes()

    err = str(capsys.readouterr().err)
    assert "unrecognized field 'y__lab'" in err
    assert (
        "'xlab': expected type 'Optional[str]', got 'bool' True" in err
        or "'xlab': expected type 'Union[str, NoneType]', got 'bool' True" in err
    )
    assert "'ymin': expected type 'Union[float, int, NoneType]', got 'str' '0'" in err

    if not strict:
        # Still should produce output unless strict mode:
        assert len(report.plot_by_id) == 1
        assert f"{id}-section-plot" in report.plot_by_id
        assert report.plot_by_id[f"{id}-section-plot"].id == f"{id}-section-plot"
        assert report.plot_by_id[f"{id}-section-plot"].plot_type == "xy_line"
        assert report.plot_by_id[f"{id}-section-plot"].pconfig.title == "DupRadar General Linear Model"
        assert report.plot_by_id[f"{id}-section-plot"].pconfig.xlog is True
        assert report.plot_by_id[f"{id}-section-plot"].pconfig.xlab is None  # wrong type
        assert report.plot_by_id[f"{id}-section-plot"].pconfig.ymax == 100
        assert report.plot_by_id[f"{id}-section-plot"].pconfig.ymin is None  # wrong type


def test_missing_id_and_title(tmp_path, capsys):
    id = "mysample"
    file = tmp_path / f"{id}_mqc.txt"
    file.write_text(
        """\
#plot_type: 'linegraph'
#pconfig:
#    xlab: "expression"
0.561167227833894 0.0146313784854042
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])

    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert f"{id}-section-plot" in report.plot_by_id
    assert report.plot_by_id[f"{id}-section-plot"].id == f"{id}-section-plot"
    assert report.plot_by_id[f"{id}-section-plot"].plot_type == "xy_line"
    assert report.plot_by_id[f"{id}-section-plot"].pconfig.xlab == "expression"


def test_with_separate_config(tmp_path, capsys):
    file = tmp_path / "mysample-concordance.txt"
    file.write_text("""Sample	'08021342'	'08027127'\n'08021342'	1.0	0.378""")

    conf_file = tmp_path / "multiqc_config.yaml"
    conf_file.write_text(
        """\
custom_data:
    concordance:
        id: 'concordance'
        section_name: 'Concordance Rates'
        plot_type: 'heatmap'
        pconfig:
            id: 'concordance_heatmap'
        sort_rows: true
sp:
    concordance:
        fn: '*concordance.txt'
"""
    )

    report.analysis_files = [file]
    update_config(cfg=ClConfig(config_files=[conf_file]))

    file_search()
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert "concordance_heatmap" in report.plot_by_id
    assert report.plot_by_id["concordance_heatmap"].id == "concordance_heatmap"
    assert report.plot_by_id["concordance_heatmap"].plot_type == "heatmap"
    assert len(report.plot_by_id["concordance_heatmap"].datasets) == 1
    assert report.plot_by_id["concordance_heatmap"].datasets[0].rows == [[1.0, 0.378]]
    assert report.plot_by_id["concordance_heatmap"].datasets[0].xcats == ["08021342", "08027127"]
    assert report.plot_by_id["concordance_heatmap"].datasets[0].ycats == ["08021342"]


def test_full_run_with_config(tmp_path, capsys):
    file = tmp_path / "mysample-concordance.txt"
    file.write_text("""Sample	'08021342'	'08027127'\n'08021342'	1.0	0.378""")

    conf_file = tmp_path / "multiqc_config.yaml"
    conf_file.write_text(
        """\
custom_data:
    concordance:
        id: 'concordance'
        section_name: 'Concordance Rates'
        plot_type: 'heatmap'
        pconfig:
            id: 'concordance_heatmap'
        sort_rows: true
sp:
    concordance:
        fn: '*concordance.txt'
"""
    )

    multiqc.run(  # pylint: disable=no-member
        file,
        cfg=ClConfig(
            strict=True,
            force=True,
            config_files=[conf_file],
            filename="stdout",
        ),
    )

    out = capsys.readouterr().out
    assert '<h2 class="mqc-module-title" id="concordance">Concordance Rates</h2>' in out
    assert '<div class="mqc-section mqc-section-concordance">' in out
    assert 'value="0.378"' in out


def test_mqc_ext_match_custom_op(tmp_path):
    """
    Verify that _mqc.tsv match the specific patterns in sp before being discovered by generic "custom_content" patterns
    """
    file1 = tmp_path / "target___test1.o2o_aln.tsv"
    file2 = tmp_path / "target___test2.o2o_aln_mqc.tsv"
    file1.write_text(
        """\
Sample	Metric
target___test1	1
"""
    )
    file2.write_text(
        """\
Sample	Metric
target___test2	2
"""
    )

    conf = tmp_path / "multiqc_config.yaml"
    conf.write_text("""
custom_data:
  last_o2o:
    plot_type: "table"

sp:
  last_o2o:
    fn: "target__*tsv"
""")

    report.analysis_files = [file1, file2]
    update_config(cfg=ClConfig(config_files=[conf], run_modules=["custom_content"]))

    file_search()
    custom_module_classes()

    # Expecting to see only one table, and no bar plot from the _mqc file
    assert len(report.plot_by_id) == 1
    assert "last_o2o-section-plot" in report.plot_by_id
    assert report.plot_by_id["last_o2o-section-plot"].id == "last_o2o-section-plot"
    assert report.plot_by_id["last_o2o-section-plot"].plot_type == "violin"


@pytest.mark.parametrize(
    ["section_name", "is_good", "contents"],
    [
        (
            None,
            True,
            """\
FILE	SEQUENCE	START	END	STRAND	GENE
myfile.fasta	chr1	55312	56664	+	GENE""",
        ),
        (
            None,
            True,
            """\
#FILE	SEQUENCE	START	END	STRAND	GENE
myfile.fasta	chr1	55312	56664	+	GENE""",
        ),
        (
            "My section",
            True,
            """\
#section_name: "My section"
#FILE	SEQUENCE	START	END	STRAND	GENE
myfile.fasta	chr1	55312	56664	+	GENE""",
        ),
        (
            "My section",
            True,
            """\
#section_name: "My section"
FILE	SEQUENCE	START	END	STRAND	GENE
myfile.fasta	chr1	55312	56664	+	GENE""",
        ),
        (
            None,
            True,
            """\
#key: value
FILE	SEQUENCE	START	END	STRAND	GENE
myfile.fasta	chr1	55312	56664	+	GENE""",
        ),
        (
            None,
            False,
            """\
#section_name: "Missing closing quote
FILE	SEQUENCE	START	END	STRAND	GENE
myfile.fasta	chr1	55312	56664	+	GENE""",
        ),
    ],
)
def test_from_tsv(tmp_path, section_name, is_good, contents):
    tmp_path.joinpath("mysample_mqc.tsv").write_text(contents)

    report.analysis_files = [tmp_path]
    update_config(cfg=ClConfig(run_modules=["custom_content"]))

    file_search()
    if not is_good:
        with pytest.raises(ConfigValidationError):
            custom_module_classes()
        return

    custom_module_classes()
    assert len(report.plot_by_id) == 1
    assert "mysample-section-plot" in report.plot_by_id
    assert report.plot_by_id["mysample-section-plot"].plot_type == "violin"
    assert len(report.plot_by_id["mysample-section-plot"].datasets) == 1
    assert report.plot_by_id["mysample-section-plot"].datasets[0].header_by_metric.keys() == {
        "SEQUENCE",
        "START",
        "END",
        "STRAND",
        "GENE",
    }

    assert report.plot_by_id["mysample-section-plot"].datasets[0].violin_value_by_sample_by_metric == {
        "SEQUENCE": {"myfile.fasta": "chr1"},
        "START": {"myfile.fasta": 55312.0},
        "END": {"myfile.fasta": 56664.0},
        "STRAND": {"myfile.fasta": "+"},
        "GENE": {"myfile.fasta": "GENE"},
    }
    assert report.plot_by_id["mysample-section-plot"].layout.title.text == "My section" if section_name else "mysample"


def test_heatmap_with_numerical_cats(tmp_path):
    plot_id = "my_plot"
    file = tmp_path / "mysample_mqc.json"
    file.write_text(
        f"""\
{{
    "id": "{plot_id}",
    "plot_type": "heatmap",
    "pconfig": {{
        "title": "Annotation stats (DRAMv)",
        "min": 0
    }},
    "ycats": ["sample 1", "sample 2"],
    "xcats": [1, 2, 3, 4],
    "data": [[0.9, 0.87, 0.73, 0], [0, 1, 0, 0.7]]
}}
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert f"{plot_id}-section-plot" in report.plot_by_id
    assert report.plot_by_id[f"{plot_id}-section-plot"].id == f"{plot_id}-section-plot"
    assert report.plot_by_id[f"{plot_id}-section-plot"].plot_type == "heatmap"


def test_on_all_example_files(data_dir):
    """
    Run on all example in data/custom_content, verify it didn't fail.
    Deprecate this in the future in favour of more granular tests like those above.
    """
    report.analysis_files = [data_dir]
    config.run_modules = ["custom_content"]

    file_search()
    custom_module_classes()
