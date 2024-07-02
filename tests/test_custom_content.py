import pytest

from multiqc import report
from multiqc.modules.custom_content import custom_module_classes


def test_custom_content(data_dir):
    subdir = data_dir / "custom_content" / "line_plot"

    report.reset()
    report.analysis_files = [subdir]
    report.search_files(["custom_content"])
    custom_module_classes()

    assert len(report.plot_data) == 1
    assert "dupradar" in report.plot_data
    assert report.plot_data["dupradar"]["id"] == "dupradar"
    assert report.plot_data["dupradar"]["plot_type"] == "linegraph"
