from typing import Any, Dict

import pytest

from multiqc import config, report


@pytest.fixture
def run_fgumi(tmp_path):
    """Factory: write ``{filename: content}`` to a temp dir and run the fgumi module on those files."""
    from multiqc.modules.fgumi import MultiqcModule

    # Module tests live outside tests/, so tests/conftest.py's autouse reset does not reach them. Reset here, before
    # the test body, so config a test sets before calling run_fgumi (sample filters, fgumi_config) is kept.
    config.reset()

    def _run(files):
        paths = []
        for name, content in files.items():
            path = tmp_path / name
            path.write_text(content)
            paths.append(path)
        report.reset()
        config.preserve_module_raw_data = True
        report.analysis_files = paths
        report.search_files(["fgumi"])
        return MultiqcModule()

    yield _run
    config.reset()


def line_points(plot_id: str, tab: int, sample: str) -> Dict[Any, Any]:
    """The ``{x: y}`` points drawn for ``sample`` on tab ``tab`` of line plot ``plot_id``."""
    lines = [line for line in report.plot_by_id[plot_id].datasets[tab].lines if line.name == sample]
    return dict(lines[0].pairs) if lines else {}


def general_stats(sample: str) -> Dict[str, Any]:
    """Every General Statistics value the module added for ``sample``, across all its column groups."""
    values: Dict[str, Any] = {}
    for rows_by_sample in report.general_stats_data.values():
        for row in rows_by_sample.get(sample, []):
            values.update(row.data)
    return values
