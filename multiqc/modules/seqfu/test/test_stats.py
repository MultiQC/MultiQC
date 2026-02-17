import pytest

from multiqc import report, config
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.types import SampleName
from multiqc.modules.seqfu.stats import all_same_length
from multiqc.modules.seqfu.seqfu import MultiqcModule
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


NUM_SAMPLES_PER_FILE = {
    "paired_end.tsv": 50,
    "promethion_fail.tsv": 1,
    "promethion_pass_stdin.tsv": 1,
    "promethion_pass.tsv": 1,
    "single_end_gc.tsv": 7,
    "single_end.tsv": 7,
}


def test_data_parsed(data_dir):
    data_subdir = data_dir / "modules/seqfu/stats"
    assert data_subdir.exists()
    for path in data_subdir.rglob("*.tsv"):
        print(path)
        path = data_subdir / path
        if path.name not in NUM_SAMPLES_PER_FILE:
            continue

        report.reset()
        report.analysis_files = [path]
        report.search_files(["seqfu"])
        config.preserve_module_raw_data = True
        m = MultiqcModule()
        assert m.saved_raw_data is not None
        assert len(list(m.saved_raw_data.values())[0]) == NUM_SAMPLES_PER_FILE[path.name]


def test_empty_file_parsing(data_dir):
    path = data_dir / "modules/seqfu/stats/empty.tsv"
    assert path.exists()
    report.reset()
    report.analysis_files = [path]
    report.search_files(["seqfu"])
    with pytest.raises(ModuleNoSamplesFound):
        _ = MultiqcModule()


def test_use_filename_as_sample_name(data_dir):
    path = data_dir / "modules/seqfu/stats/promethion_pass.tsv"
    assert path.exists()
    report.reset()
    report.analysis_files = [path]
    report.search_files(["seqfu"])
    config.preserve_module_raw_data = True
    config.use_filename_as_sample_name = True
    m = MultiqcModule()
    assert m.saved_raw_data is not None
    assert len(list(m.saved_raw_data.values())[0]) == NUM_SAMPLES_PER_FILE[path.name]
    assert m._clean_s_name(path.name) in list(m.saved_raw_data.values())[0]


def test_all_same_length():
    assert all_same_length({SampleName("a"): {"Min": 1, "Max": 1}, SampleName("b"): {"Min": 1, "Max": 1}})
    assert not all_same_length({SampleName("a"): {"Min": 1, "Max": 1}, SampleName("b"): {"Min": 1, "Max": 2}})
