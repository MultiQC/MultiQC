import pytest
from multiqc.utils import testing
from multiqc.utils.testing import BaseModuleTest
from multiqc.modules.fastqc.fastqc import MultiqcModule

@pytest.fixture
def data_dir():
    return testing.data_dir()

class TestFastqcModule(BaseModuleTest):
    MODULE_CLASS = MultiqcModule
    MODULE_NAME = "fastqc"

    @pytest.mark.parametrize("subdir", [
        "",  # root fastqc_data.txt
        "v0.11.2",  # zipped report
        "zero_reads",  # zero reads
        "nan_reads",  # nan reads
        "groups",  # group/multiple samples
        "issue_1941_unicode",  # unicode sample names
        "issue_2343_numerical_samples",  # numerical sample names
    ])
    def test_fastqc_snapshot(self, data_dir, snapshot, subdir):
        module_data_dir = data_dir / "modules" / self.MODULE_NAME / subdir
        # Find all fastqc_data.txt and *_fastqc.zip in the directory
        data_files = list(module_data_dir.glob("fastqc_data.txt")) + list(module_data_dir.glob("*_fastqc.zip"))
        if not data_files:
            pytest.skip(f"No fastqc test files found in {module_data_dir}")
        module_snapshot = self.create_module_snapshot(data_dir, snapshot_config={})
        self.assert_module_data_integrity(module_snapshot)
        snapshot.assert_match(module_snapshot.get_complete_snapshot()) 