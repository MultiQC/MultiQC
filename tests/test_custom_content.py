import tempfile

import pytest

import multiqc
from multiqc import config, report
from multiqc.core.file_search import file_search
from multiqc.core.update_config import ClConfig, update_config
from multiqc.modules.custom_content import custom_module_classes
from multiqc.types import Anchor, ColumnKey, SampleGroup
from multiqc.validation import ModuleConfigValidationError


def test_linegraph_single_sample_txt(data_dir):
    path = data_dir / "custom_content" / "embedded_config" / "linegraph_single_sample_txt_mqc.txt"
    """
#id: dupradar
#plot_type: 'linegraph'
#section_name: 'DupRadar'
#section_href: 'bioconductor.org/packages/release/bioc/html/dupRadar.html'
#description: "provides duplication rate quality control for RNA-Seq datasets. Highly expressed genes can be expected to have a lot of duplicate reads, but high numbers of duplicates at low read counts can indicate low library complexity with technical duplication.
#    This plot shows the general linear models - a summary of the gene duplication distributions. "
#pconfig:
#    title: 'DupRadar General Linear Model'
#    xlog: True
#    xlab: 'expression (reads/kbp)'
#    ylab: '% duplicate reads'
#    ymax: 100
#    ymin: 0
#    tt_label: '<b>{point.x:.1f} reads/kbp</b>: {point.y:,.2f}% duplicates'
#    x_lines:
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

    report.analysis_files = [path]
    report.search_files(["custom_content"])
    custom_module_classes()

    anchor = Anchor("dupradar-section-plot")
    assert len(report.plot_by_id) == 1
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == "dupradar"
    assert report.plot_by_id[anchor].plot_type == "xy_line"


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
    anchor = Anchor(f"{id}-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == id
    assert report.plot_by_id[anchor].plot_type == "xy_line"

    err = str(capsys.readouterr().err)
    assert "Line plot's x_lines or y_lines 'label' field is expected to be a string" in err
    assert "'LongDash' is a deprecated dash style, use 'longdash'" in err
    assert "deprecated field 'colour'. Use 'color' instead" in err
    assert "deprecated field 'xLog'. Use 'xlog' instead" in err
    assert "deprecated field 'xPlotLines'. Use 'x_lines' instead" in err


@pytest.mark.parametrize("strict", [True, False])
def test_wrong_fields(tmp_path, caplog, strict, monkeypatch):
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
#    y__lab: '% duplicate reads'
#    xlab: True
#    xlog: 'True'
#    ymax: [1, 2]
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
        with pytest.raises(ModuleConfigValidationError):
            custom_module_classes()
    else:
        custom_module_classes()

    out = caplog.text
    assert "unrecognized field 'y__lab'" in out
    assert (
        "• 'ymax': expected type '<class 'float'>', got 'list' [1, 2]" in out
        or "• 'ymax': expected type 'Union[float, int, NoneType]', got 'list' [1, 2]" in out  # Python < 3.10
    )

    if not strict:
        # Still should produce output unless strict mode:
        assert len(report.plot_by_id) == 1
        anchor = Anchor(f"{id}-section-plot")
        assert anchor in report.plot_by_id
        assert report.plot_by_id[anchor].id == id
        assert report.plot_by_id[anchor].plot_type == "xy_line"
        assert report.plot_by_id[anchor].pconfig.title == "DupRadar General Linear Model"
        assert (
            report.plot_by_id[anchor].pconfig.xlab == "True"  # cast to string
            or report.plot_by_id[anchor].pconfig.xlab is None  # Python < 3.10
        )
        assert report.plot_by_id[anchor].pconfig.xlog is True
        assert report.plot_by_id[anchor].pconfig.ymax is None
        assert (
            report.plot_by_id[anchor].pconfig.ymin == 0  # cast to int
            or report.plot_by_id[anchor].pconfig.ymin is None  # Python < 3.10
        )


def test_missing_id_and_title(tmp_path):
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
    anchor = Anchor(f"{id}-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == id
    assert report.plot_by_id[anchor].plot_type == "xy_line"
    assert report.plot_by_id[anchor].pconfig.xlab == "expression"


def test_with_separate_config_and_quotes(data_dir):
    file = data_dir / "custom_content" / "with_config" / "run_concordance" / "run_concordance.txt"
    """
Sample	'08021342'	'08027127'
'08021342'	1.0	0.378
'08027127'	0.341	1.0
    """
    conf_file = data_dir / "custom_content" / "with_config" / "run_concordance" / "multiqc_config.yaml"
    """
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

    report.analysis_files = [file]
    update_config(cfg=ClConfig(config_files=[conf_file]))

    file_search()
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    anchor = Anchor("concordance_heatmap")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == anchor
    assert report.plot_by_id[anchor].plot_type == "heatmap"
    assert len(report.plot_by_id[anchor].datasets) == 1
    assert report.plot_by_id[anchor].datasets[0].rows[0][0] == 1.0
    assert report.plot_by_id[anchor].datasets[0].rows[0][1] == 0.378
    assert report.plot_by_id[anchor].datasets[0].xcats[0] == "08021342"
    assert report.plot_by_id[anchor].datasets[0].ycats[0] == "08021342"


def test_full_run_with_config(data_dir, capsys):
    file = data_dir / "custom_content" / "with_config" / "run_concordance" / "run_concordance.txt"
    conf_file = data_dir / "custom_content" / "with_config" / "run_concordance" / "multiqc_config.yaml"

    multiqc.run(  # pylint: disable=no-member
        file,
        cfg=ClConfig(
            strict=True,
            force=True,
            config_files=[conf_file],
            filename="stdout",
            development=True,
        ),
    )

    out = capsys.readouterr().out
    assert '<h2 class="mqc-module-title" id="concordance">Concordance Rates</h2>' in out
    assert '<div class="mqc-section mqc-section-concordance">' in out

    assert len(report.plot_by_id) == 1
    anchor = Anchor("concordance_heatmap")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == "concordance_heatmap"
    assert report.plot_by_id[anchor].plot_type == "heatmap"
    assert len(report.plot_by_id[anchor].datasets) == 1
    assert report.plot_by_id[anchor].datasets[0].rows[0][0] == 1.0
    assert report.plot_by_id[anchor].datasets[0].rows[0][1] == 0.378
    assert report.plot_by_id[anchor].datasets[0].xcats[0] == "08021342"
    assert report.plot_by_id[anchor].datasets[0].ycats[0] == "08021342"


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
    conf.write_text(
        """
        custom_data:
          last_o2o:
            plot_type: "table"
        
        sp:
          last_o2o:
            fn: "target__*tsv"
        """
    )

    report.analysis_files = [file1, file2]
    update_config(cfg=ClConfig(config_files=[conf], run_modules=["custom_content"]))

    file_search()
    custom_module_classes()

    # Expecting to see only one table, and no bar plot from the _mqc file
    assert len(report.plot_by_id) == 1
    anchor = Anchor("last_o2o-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == "last_o2o"
    assert report.plot_by_id[anchor].plot_type == "violin"


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
    id = "mysample"
    tmp_path.joinpath(f"{id}_mqc.tsv").write_text(contents)

    report.analysis_files = [tmp_path]
    update_config(cfg=ClConfig(run_modules=["custom_content"]))

    file_search()
    if not is_good:
        with pytest.raises(ModuleConfigValidationError):
            custom_module_classes()
        return

    custom_module_classes()
    assert len(report.plot_by_id) == 1
    anchor = Anchor(f"{id}-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].plot_type == "violin"
    assert len(report.plot_by_id[anchor].datasets) == 1
    assert report.plot_by_id[anchor].datasets[0].header_by_metric.keys() == {
        "SEQUENCE",
        "START",
        "END",
        "STRAND",
        "GENE",
    }

    assert report.plot_by_id[anchor].datasets[0].violin_value_by_sample_by_metric == {
        "SEQUENCE": {"myfile.fasta": "chr1"},
        "START": {"myfile.fasta": 55312.0},
        "END": {"myfile.fasta": 56664.0},
        "STRAND": {"myfile.fasta": "+"},
        "GENE": {"myfile.fasta": "GENE"},
    }
    if section_name:
        assert report.plot_by_id[anchor].layout.title.text == section_name
    else:
        assert report.plot_by_id[anchor].layout.title.text == id.title()


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
    anchor = Anchor(f"{plot_id}-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == plot_id
    assert report.plot_by_id[anchor].plot_type == "heatmap"

    assert report.plot_by_id[anchor].datasets[0].rows == [[0.9, 0.87, 0.73, 0], [0, 1, 0, 0.7]]
    assert report.plot_by_id[anchor].datasets[0].xcats == ["1", "2", "3", "4"]
    assert report.plot_by_id[anchor].datasets[0].ycats == ["sample 1", "sample 2"]


def test_on_all_example_files(data_dir):
    """
    Run on all example in data/custom_content, verify it didn't fail.
    Deprecate this in the future in favour of more granular tests like those above.
    """

    report.analysis_files = [data_dir / "custom_content"]
    config.run_modules = ["custom_content"]

    file_search()
    custom_module_classes()


def test_boxplot_json(tmp_path):
    file = tmp_path / "mysample_mqc.json"
    file.write_text(
        """\
{
    "id": "boxplot",
    "plot_type": "box",
    "pconfig": {
        "title": "Boxplot"
    },
    "data": {
        "sample 1": [1, 2, 3, 4],
        "sample 2": [1.1, 0.2, 3.3, 4.4, 5.5]
    }
}
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    anchor = Anchor("boxplot-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].id == "boxplot"
    assert report.plot_by_id[anchor].plot_type == "box"
    assert len(report.plot_by_id[anchor].datasets) == 1
    assert report.plot_by_id[anchor].datasets[0].data[
        report.plot_by_id[anchor].datasets[0].samples.index("sample 1")
    ] == [1, 2, 3, 4]
    assert report.plot_by_id[anchor].datasets[0].data[
        report.plot_by_id[anchor].datasets[0].samples.index("sample 2")
    ] == [1.1, 0.2, 3.3, 4.4, 5.5]


def test_boxplot_txt(tmp_path):
    file = tmp_path / "mysample_mqc.txt"
    file.write_text(
        """\
#plot_type: 'box'
Sample1	Sample2	Sample3
1.2	2.3	3.4
2.1	3.2	4.3
1.5	2.5	3.5
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    anchor = Anchor("mysample-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].plot_type == "box"
    assert len(report.plot_by_id[anchor].datasets[0].data) == 3  # 3 groups
    assert report.plot_by_id[anchor].datasets[0].data == [[3.4, 4.3, 3.5], [2.3, 3.2, 2.5], [1.2, 2.1, 1.5]]


def test_scatter_plot_parsing(tmp_path):
    """Test scatter plot data parsing"""
    file = tmp_path / "scatter_mqc.txt"
    file.write_text(
        """\
#plot_type: 'scatter'
A	1.2	2.3
B	2.1	3.2
C	1.5	2.5
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    anchor = Anchor("scatter-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].plot_type == "scatter"
    assert len(report.plot_by_id[anchor].datasets[0].points) == 3  # 3 samples
    assert report.plot_by_id[anchor].datasets[0].points[0] == {"x": 1.2, "y": 2.3, "name": "A"}


def test_line_plot_with_x_labels(tmp_path):
    """Test line plot with x-axis labels in first row"""
    file = tmp_path / "line_mqc.txt"
    file.write_text(
        """\
#plot_type: 'linegraph'
	t0	t1	t2	t3
sample1	1.0	1.2	1.1	0.9
sample2	0.8	1.4	1.6	1.2
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    anchor = Anchor("line-section-plot")
    assert anchor in report.plot_by_id
    assert report.plot_by_id[anchor].plot_type == "xy_line"
    assert len(report.plot_by_id[anchor].datasets[0].lines) == 2  # 2 samples
    assert report.plot_by_id[anchor].datasets[0].lines[0].name == "sample1"
    assert report.plot_by_id[anchor].datasets[0].lines[0].pairs == [("t0", 1.0), ("t1", 1.2), ("t2", 1.1), ("t3", 0.9)]


def test_general_stats_with_config(tmp_path):
    """Test general stats with custom headers config"""
    file = tmp_path / "stats_mqc.txt"
    file.write_text(
        """\
#plot_type: 'generalstats'
#pconfig:
#    namespace: 'My Stats'
#headers:
#    value1:
#        title: 'Value 1'
#        description: 'First value'
#        max: 100
#    value2:
#        title: 'Value 2'
#        description: 'Second value'
#        min: 0
Sample	value1	value2
A	85	42
B	90	38
C	95	45
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    # General stats are added directly to the report
    assert len(report.general_stats_data) > 0
    assert "A" in report.general_stats_data[-1]
    assert report.general_stats_data[-1][SampleGroup("A")][0].data[ColumnKey("value1")] == 85
    assert report.general_stats_data[-1][SampleGroup("A")][0].data[ColumnKey("value2")] == 42
    assert report.general_stats_headers[-1][ColumnKey("value1")]["title"] == "Value 1"
    assert report.general_stats_headers[-1][ColumnKey("value2")]["title"] == "Value 2"


def test_quoted_strings_handling(tmp_path):
    """Test handling of quoted strings in data"""
    file = tmp_path / "quoted_mqc.txt"
    file.write_text(
        """\
#plot_type: 'table'
Sample	"Group Name"	'Value'
"Sample 1"	"Group A"	4_2
'Sample 2'	'Group B'	"3_8"
Sample 3	Group C	'4_5'
"""
    )

    report.analysis_files = [file]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_by_id) == 1
    anchor = Anchor("quoted-section-plot")
    assert anchor in report.plot_by_id
    plot = report.plot_by_id[anchor]
    assert plot.plot_type == "violin"

    # Check that quotes were properly stripped
    assert plot.datasets[0].all_samples == ["Sample 1", "Sample 2", "Sample 3"]
    assert plot.datasets[0].header_by_metric.keys() == {"Group_Name", "Value"}
    assert plot.datasets[0].violin_value_by_sample_by_metric == {
        "Group_Name": {"Sample 1": "Group A", "Sample 2": "Group B", "Sample 3": "Group C"},
        "Value": {"Sample 1": 42.0, "Sample 2": "3_8", "Sample 3": "4_5"},
    }
    assert plot.datasets[0].dt.sections[0].column_by_key.keys() == {"Group Name", "Value"}
    assert plot.datasets[0].dt.sections[0].column_by_key["Group Name"].title == "Group Name"
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup.keys() == {"Sample 1", "Sample 2", "Sample 3"}
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup["Sample 1"][0].sample == "Sample 1"
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup["Sample 1"][0].raw_data == {
        "Group Name": "Group A",
        "Value": 42.0,
    }
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup["Sample 2"][0].sample == "Sample 2"
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup["Sample 2"][0].raw_data == {
        "Group Name": "Group B",
        "Value": "3_8",
    }
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup["Sample 3"][0].sample == "Sample 3"
    assert plot.datasets[0].dt.sections[0].rows_by_sgroup["Sample 3"][0].raw_data == {
        "Group Name": "Group C",
        "Value": "4_5",
    }
