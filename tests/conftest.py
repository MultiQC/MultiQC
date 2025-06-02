import pytest

from multiqc import report, config, validation
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


@pytest.fixture(autouse=True)
def reset():
    """
    Reset MultiQC session after use: reset config and report
    """

    yield

    report.reset()
    config.reset()
    validation.reset()


@pytest.fixture
def multiqc_reset():
    """
    Explicit reset fixture for tests that need to reset state mid-test.
    """
    report.reset()
    config.reset()
    validation.reset()


@pytest.fixture
def snapshot_config():
    """
    Configuration for snapshot testing.

    Returns a dictionary with common configuration options for snapshot tests.
    """
    return {
        "preserve_module_raw_data": True,
        "strict": True,
        "make_data_dir": False,  # Don't create data directories during testing
    }


class BaseModuleTest:
    """
    Base class for module snapshot tests.

    This class provides common functionality for testing MultiQC modules
    with snapshot testing capabilities.
    """

    # Override these in subclasses
    MODULE_CLASS = None
    MODULE_NAME = None
    DATA_DIR_NAME = None  # If different from MODULE_NAME

    @pytest.fixture
    def module_data_dir(self, data_dir):
        """Get the data directory for this module."""
        dir_name = self.DATA_DIR_NAME or self.MODULE_NAME
        return data_dir / "modules" / dir_name

    @pytest.fixture
    def module_snapshot(self, module_data_dir, snapshot_config):
        """
        Create a snapshot of the module with all available test data.

        This fixture automatically finds all test files in the module's data directory
        and runs the module on them, returning a ModuleSnapshot instance.
        """
        if self.MODULE_CLASS is None:
            pytest.skip("MODULE_CLASS not defined")

        # Apply snapshot configuration
        for key, value in snapshot_config.items():
            setattr(config, key, value)

        # Run the module test
        return testing.run_module_test(
            module_class=self.MODULE_CLASS,
            data_files=list(module_data_dir.rglob("*")),
            config_updates=snapshot_config,
        )

    def test_module_data_integrity(self, module_snapshot):
        """Test basic data integrity for the module."""
        testing.assert_module_data_integrity(module_snapshot)

    def test_module_complete_snapshot(self, module_snapshot, snapshot):
        """Test the complete module snapshot."""
        complete_snapshot = module_snapshot.get_complete_snapshot()
        assert complete_snapshot == snapshot

    def test_module_raw_data_snapshot(self, module_snapshot, snapshot):
        """Test just the raw data snapshot."""
        raw_data = module_snapshot.get_saved_raw_data()
        assert raw_data == snapshot

    def test_module_general_stats_snapshot(self, module_snapshot, snapshot):
        """Test the general statistics snapshot."""
        general_stats = {
            "data": module_snapshot.get_general_stats_data(),
            "headers": module_snapshot.get_general_stats_headers(),
        }
        assert general_stats == snapshot

    def test_module_sections_snapshot(self, module_snapshot, snapshot):
        """Test the sections snapshot."""
        sections = module_snapshot.get_sections_data()
        assert sections == snapshot


def pytest_configure(config):
    """Configure pytest with custom markers for snapshot testing."""
    config.addinivalue_line("markers", "snapshot: mark test as a snapshot test")
    config.addinivalue_line("markers", "module_test: mark test as a module test")
