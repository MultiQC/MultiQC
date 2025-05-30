import pytest
from multiqc.utils import testing
from multiqc.utils.testing import BaseModuleTest
from multiqc.modules.fastqc.fastqc import MultiqcModule

@pytest.fixture
def data_dir():
    return testing.data_dir()

@pytest.fixture
def snapshot_config():
    """Configuration for snapshot testing."""
    return {
        "preserve_module_raw_data": True,
        "strict": True,
        "make_data_dir": False,
    }

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
    def test_fastqc_snapshot(self, data_dir, snapshot, snapshot_config, subdir):
        module_data_dir = data_dir / "modules" / self.MODULE_NAME / subdir
        # Find all fastqc_data.txt and *_fastqc.zip in the directory
        data_files = list(module_data_dir.glob("fastqc_data.txt")) + list(module_data_dir.glob("*_fastqc.zip"))
        if not data_files:
            pytest.skip(f"No fastqc test files found in {module_data_dir}")
        
        # Run the module test on the specific subdirectory
        module_snapshot = testing.run_module_test(
            module_class=self.MODULE_CLASS,
            data_files=data_files,
            config_updates=snapshot_config,
        )
        
        # Assert basic data integrity
        self.assert_module_data_integrity(module_snapshot)
        
        # Snapshot the raw data output from the module
        raw_data = module_snapshot.get_saved_raw_data()
        assert raw_data == snapshot 