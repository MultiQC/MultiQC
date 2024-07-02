import pytest

from multiqc import report
from multiqc.core.update_config import update_config, ClConfig
from multiqc.modules.custom_content import custom_module_classes
from multiqc.validation import ConfigValidationError


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

    report.reset()
    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert f"{id}-plot" in report.plot_by_id
    assert report.plot_by_id[f"{id}-plot"].id == f"{id}-plot"
    assert report.plot_by_id[f"{id}-plot"].plot_type == "xy_line"


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

    report.reset()
    report.analysis_files = [file]
    update_config(cfg=ClConfig(verbose=True))
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    assert f"{id}-plot" in report.plot_by_id
    assert report.plot_by_id[f"{id}-plot"].id == f"{id}-plot"
    assert report.plot_by_id[f"{id}-plot"].plot_type == "xy_line"

    err = str(capsys.readouterr().err)
    assert "Line plot's x_lines or y_lines 'label' field is expected to be a string" in err
    assert "Deprecated field 'xLog'. Use 'xlog' instead" in err
    assert "Deprecated field 'xPlotLines'. Use 'x_lines' instead" in err


@pytest.mark.parametrize("strict", [True, False])
def test_wrong_fields(tmp_path, capsys, strict):
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

    report.reset()
    report.analysis_files = [file]
    report.search_files(["custom_content"])
    update_config(cfg=ClConfig(strict=strict))

    if strict:
        with pytest.raises(ConfigValidationError):
            custom_module_classes()
    else:
        custom_module_classes()

    err = str(capsys.readouterr().err)
    assert "unrecognized field 'y__lab'" in err
    assert "'xlab': expected type 'Optional[str]', got 'bool' True" in err
    assert "'ymin': expected type 'Union[float, int, NoneType]', got 'str' '0'" in err

    if not strict:
        # Still should produce output unless strict mode:
        assert len(report.plot_by_id) == 1
        assert f"{id}-plot" in report.plot_by_id
        assert report.plot_by_id[f"{id}-plot"].id == f"{id}-plot"
        assert report.plot_by_id[f"{id}-plot"].plot_type == "xy_line"
        assert report.plot_by_id[f"{id}-plot"].pconfig.title == "DupRadar General Linear Model"
        assert report.plot_by_id[f"{id}-plot"].pconfig.xlog is True
        assert report.plot_by_id[f"{id}-plot"].pconfig.xlab is None  # wrong type
        assert report.plot_by_id[f"{id}-plot"].pconfig.ymax == 100
        assert report.plot_by_id[f"{id}-plot"].pconfig.ymin is None  # wrong type


@pytest.mark.parametrize("strict", [True, False])
def test_missing_id_and_title(tmp_path, capsys, strict):
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

    report.reset()
    report.analysis_files = [file]
    report.search_files(["custom_content"])
    update_config(cfg=ClConfig(strict=strict))

    custom_module_classes()
    err = str(capsys.readouterr().err)
    print(err)

    # Still should produce output unless strict mode:
    assert len(report.plot_by_id) == 1
    assert f"{id}-plot" in report.plot_by_id
    assert report.plot_by_id[f"{id}-plot"].id == f"{id}-plot"
    assert report.plot_by_id[f"{id}-plot"].plot_type == "xy_line"
    assert report.plot_by_id[f"{id}-plot"].pconfig.xlab == "expression"
