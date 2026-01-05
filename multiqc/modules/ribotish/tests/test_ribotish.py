"""Tests for the Ribo-TISH module."""

import pytest

from multiqc import report
from multiqc.base_module import ModuleNoSamplesFound
from multiqc.modules.ribotish import MultiqcModule


# Valid qual file content - line 4 has the frame counts dict
VALID_QUAL_BASIC = """\
# comment line 1
# comment line 2
# comment line 3
{28: [100, 20, 30], 29: [150, 25, 25], 30: [80, 10, 10]}
"""

VALID_QUAL_WITH_DOLLAR = """\
# comment line 1
# comment line 2
# comment line 3
$ {28: [100, 20, 30], 29: [150, 25, 25]}
"""

VALID_QUAL_SINGLE_LENGTH = """\
line1
line2
line3
{28: [100, 20, 30]}
"""

VALID_QUAL_ALL_ZEROS = """\
line1
line2
line3
{25: [0, 0, 0], 26: [0, 0, 0], 27: [0, 0, 0]}
"""

# Invalid qual file contents
INVALID_TOO_FEW_LINES = """line1
line2
line3"""

INVALID_NOT_A_DICT = """\
line1
line2
line3
this is not a dict
"""

INVALID_SYNTAX_ERROR = """\
line1
line2
line3
{25: [1, 2, 3"""

INVALID_LIST_NOT_DICT = """\
line1
line2
line3
[1, 2, 3]
"""

INVALID_VALUES_NOT_LISTS = """\
line1
line2
line3
{25: 'not a list', 26: 123}
"""

INVALID_WRONG_LIST_LENGTH = """\
line1
line2
line3
{25: [1, 2], 26: [1, 2, 3, 4]}
"""

INVALID_NON_NUMERIC = """\
line1
line2
line3
{25: ['a', 'b', 'c']}
"""


@pytest.fixture
def run_ribotish_module(tmp_path):
    """Factory to run ribotish module with custom file content."""

    def _run_module(file_content: str, filename: str = "test_qual.txt"):
        # Write file to temp path
        test_file = tmp_path / filename
        test_file.write_text(file_content)

        # Set up report to search this file
        report.reset()
        report.analysis_files = [test_file]
        report.search_files(["ribotish"])

        # Run the module
        return MultiqcModule()

    return _run_module


class TestRibotishParsing:
    """Test parsing of qual files."""

    def test_valid_basic_parsing(self, run_ribotish_module):
        """Test parsing a valid qual file."""
        module = run_ribotish_module(VALID_QUAL_BASIC)

        assert "test" in module.ribotish_data
        data = module.ribotish_data["test"]

        assert 28 in data
        assert 29 in data
        assert 30 in data
        assert data[28] == [100, 20, 30]
        assert data[29] == [150, 25, 25]
        assert data[30] == [80, 10, 10]

    def test_dollar_prefix_stripped(self, run_ribotish_module):
        """Test that $ prefix on line 4 is properly stripped."""
        module = run_ribotish_module(VALID_QUAL_WITH_DOLLAR)

        assert "test" in module.ribotish_data
        assert 28 in module.ribotish_data["test"]
        assert 29 in module.ribotish_data["test"]

    def test_single_read_length(self, run_ribotish_module):
        """Test handling of single read length."""
        module = run_ribotish_module(VALID_QUAL_SINGLE_LENGTH)

        assert "test" in module.ribotish_data
        assert len(module.ribotish_data["test"]) == 1
        assert 28 in module.ribotish_data["test"]


class TestRibotishMalformedFiles:
    """Test handling of malformed qual files."""

    def test_file_with_fewer_than_4_lines(self, run_ribotish_module):
        """Test file with <4 lines is skipped."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_TOO_FEW_LINES)

    def test_file_with_invalid_dict_format(self, run_ribotish_module):
        """Test file with invalid dict format on line 4."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_NOT_A_DICT)

    def test_file_with_syntax_error(self, run_ribotish_module):
        """Test file with Python syntax error on line 4."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_SYNTAX_ERROR)

    def test_file_with_list_instead_of_dict(self, run_ribotish_module):
        """Test file where line 4 is a list instead of dict."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_LIST_NOT_DICT)

    def test_file_with_wrong_value_type(self, run_ribotish_module):
        """Test file where dict values are not lists."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_VALUES_NOT_LISTS)

    def test_file_with_wrong_list_length(self, run_ribotish_module):
        """Test file where frame count lists don't have exactly 3 elements."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_WRONG_LIST_LENGTH)

    def test_file_with_non_numeric_counts(self, run_ribotish_module):
        """Test file where frame counts contain non-numeric values."""
        with pytest.raises(ModuleNoSamplesFound):
            run_ribotish_module(INVALID_NON_NUMERIC)


class TestRibotishFrameCalculations:
    """Test frame proportion calculations."""

    def test_frame_proportions_sum_to_one(self, run_ribotish_module):
        """Test that frame proportions sum to 1.0 for each length."""
        module = run_ribotish_module(VALID_QUAL_BASIC)

        sample_name = "test"
        for length, props in module.frame_proportions[sample_name].items():
            total_prop = props["f0_prop"] + props["f1_prop"] + props["f2_prop"]
            assert abs(total_prop - 1.0) < 0.001, f"Proportions don't sum to 1.0 for length {length}"

    def test_frame_proportions_correct_values(self, run_ribotish_module):
        """Test that frame proportions are calculated correctly."""
        module = run_ribotish_module(VALID_QUAL_SINGLE_LENGTH)

        sample_name = "test"
        props = module.frame_proportions[sample_name][28]

        # counts are [100, 20, 30], total = 150
        assert abs(props["f0_prop"] - 100 / 150) < 0.001
        assert abs(props["f1_prop"] - 20 / 150) < 0.001
        assert abs(props["f2_prop"] - 30 / 150) < 0.001
        assert props["total"] == 150

    def test_all_zeros_handled(self, run_ribotish_module):
        """Test handling of all zeros in counts."""
        module = run_ribotish_module(VALID_QUAL_ALL_ZEROS)

        sample_name = "test"
        for length in [25, 26, 27]:
            props = module.frame_proportions[sample_name][length]
            assert props["f0_prop"] == 0
            assert props["f1_prop"] == 0
            assert props["f2_prop"] == 0
            assert props["total"] == 0
