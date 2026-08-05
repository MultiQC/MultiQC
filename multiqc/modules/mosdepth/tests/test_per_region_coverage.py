import pytest

from multiqc.modules.mosdepth.mosdepth import (
    build_per_region_rows,
    parse_regions_bed_lines,
    parse_thresholds_bed_lines,
)

REGIONS_WITH_NAMES = """1\t2488047\t2488227\tTNFRSF14\t45.83
1\t2489098\t2489338\tNPHP4\t50.91
"""

REGIONS_WITHOUT_NAMES = """CM000663.2\t0\t10\t12.50
CM000663.2\t10\t20\t0.00
"""

THRESHOLDS_FULL = """#chrom\tstart\tend\tregion\t1X\t10X\t20X\t30X
1\t2488047\t2488227\tTNFRSF14\t180\t180\t180\t180
1\t2489098\t2489338\tNPHP4\t240\t240\t240\t120
"""

THRESHOLDS_EMPTY = """#chrom\tstart\tend\tregion\t20X\t30X
"""

THRESHOLDS_INVALID_HEADER = """chrom\tstart\tend\tname\t20X
1\t2488047\t2488227\tTNFRSF14\t180
"""


class TestParseRegionsBedLines:
    def test_with_name_column(self):
        result = parse_regions_bed_lines(REGIONS_WITH_NAMES.splitlines())
        assert result == {
            ("1", 2488047, 2488227): 45.83,
            ("1", 2489098, 2489338): 50.91,
        }

    def test_without_name_column(self):
        result = parse_regions_bed_lines(REGIONS_WITHOUT_NAMES.splitlines())
        assert result == {
            ("CM000663.2", 0, 10): 12.50,
            ("CM000663.2", 10, 20): 0.00,
        }

    def test_empty_file(self):
        assert parse_regions_bed_lines([]) == {}

    def test_invalid_column_count(self):
        with pytest.raises(ValueError):
            parse_regions_bed_lines(["1\t2488047\t2488227\n"])


class TestParseThresholdsBedLines:
    def test_full_output(self):
        thresholds, by_region = parse_thresholds_bed_lines(THRESHOLDS_FULL.splitlines())
        assert thresholds == [1, 10, 20, 30]
        assert by_region == {
            ("1", 2488047, 2488227): ("TNFRSF14", [180, 180, 180, 180]),
            ("1", 2489098, 2489338): ("NPHP4", [240, 240, 240, 120]),
        }

    def test_custom_threshold_levels(self):
        thresholds, by_region = parse_thresholds_bed_lines(THRESHOLDS_EMPTY.splitlines())
        assert thresholds == [20, 30]
        assert by_region == {}

    def test_invalid_header(self):
        with pytest.raises(ValueError):
            parse_thresholds_bed_lines(THRESHOLDS_INVALID_HEADER.splitlines())

    def test_empty_file_raises_value_error(self):
        with pytest.raises(ValueError):
            parse_thresholds_bed_lines([])


class TestBuildPerRegionRows:
    def test_joins_mean_coverage_by_position(self):
        thresholds = [20, 30]
        by_region = {
            ("1", 2488047, 2488227): ("TNFRSF14", [180, 180]),
        }
        mean_cov_by_region = {("1", 2488047, 2488227): 45.83}

        rows = build_per_region_rows(thresholds, by_region, mean_cov_by_region)

        assert rows == {
            "TNFRSF14": {
                "chrom": "1",
                "start": 2488047,
                "end": 2488227,
                "mean_coverage": 45.83,
                "pct_at_20x": 100.0,
                "pct_at_30x": 100.0,
            }
        }

    def test_missing_mean_coverage_omits_column(self):
        thresholds = [20]
        by_region = {("1", 0, 100): ("GENE1", [50])}

        rows = build_per_region_rows(thresholds, by_region, mean_cov_by_region={})

        assert "mean_coverage" not in rows["GENE1"]
        assert rows["GENE1"]["pct_at_20x"] == 50.0

    def test_duplicate_names_disambiguated_with_coordinates(self):
        thresholds = [20]
        by_region = {
            ("1", 0, 100): ("LRP1B", [50]),
            ("1", 200, 300): ("LRP1B", [100]),
            ("1", 400, 500): ("ATM", [100]),
        }

        rows = build_per_region_rows(thresholds, by_region, mean_cov_by_region={})

        assert set(rows.keys()) == {
            "LRP1B (1:0-100)",
            "LRP1B (1:200-300)",
            "ATM",
        }

    def test_zero_length_region_reports_zero_percent(self):
        thresholds = [20]
        by_region = {("1", 100, 100): ("GENE1", [0])}

        rows = build_per_region_rows(thresholds, by_region, mean_cov_by_region={})

        assert rows["GENE1"]["pct_at_20x"] == 0.0
