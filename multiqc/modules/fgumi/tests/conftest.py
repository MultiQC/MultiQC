import pytest

from multiqc import config, report


@pytest.fixture
def run_fgumi(tmp_path):
    """Factory: write ``{filename: content}`` to a temp dir and run the fgumi module on those files."""
    from multiqc.modules.fgumi import MultiqcModule

    original = config.preserve_module_raw_data

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
    config.preserve_module_raw_data = original
