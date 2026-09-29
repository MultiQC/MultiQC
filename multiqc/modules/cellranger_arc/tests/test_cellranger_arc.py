import glob

import pytest
from multiqc import config, report
from multiqc.modules.cellranger_arc import MultiqcModule
from multiqc.utils import testing


@pytest.fixture
def data_dir():
    return testing.data_dir()


def test_cellranger_arc(data_dir):
    html_files = glob.glob(str(data_dir / "modules/cellranger_arc/arc-*/*.html"))
    for file in html_files:
        report.analysis_files = [file]
        report.search_files(["cellranger_arc"])

    config.preserve_module_raw_data = True
    m = MultiqcModule()
    assert m.saved_raw_data is not None
    assert len(m.saved_raw_data) > 0

    ## check if both sameples are present in the data
    data = m.saved_raw_data["multiqc_cellranger_arc"]
    assert "TEST_SAMPLE" in data, "Missing sample: 'TEST_SAMPLE'"
    assert "P210059__KIZH_03" in data, "Missing sample: 'P210059__KIZH_03'"

    ## check if general stats data is there and there are 8 sections
    assert len(report.general_stats_data) > 0
    assert len(m.sections) == 8

    ## check of the plots have data
    for sam in ["TEST_SAMPLE", "P210059__KIZH_03"]:
        for key in m.plots_data_by_sample:
            assert m.plots_data_by_sample[key][sam], f"{key} plot data is empty for {sam}"


def test_atac_gex_metric_separation(data_dir):
    """
    Test that ATAC and GEX metrics are properly separated using namespace prefixing.

    This test validates the fix for metric collision between ATAC and GEX metrics that
    share the same base name (e.g., "Sequenced read pairs", "Valid barcodes").

    The fix applies namespace prefixes:
    - ATAC metrics: "ATAC_Sequenced read pairs", "ATAC_Valid barcodes", etc.
    - GEX metrics: "GEX_Sequenced read pairs", "GEX_Valid barcodes", etc.

    This ensures that ATAC and GEX values don't overwrite each other.
    """
    # Reset report state for clean test
    from multiqc import reset

    reset()

    html_files = glob.glob(str(data_dir / "modules/cellranger_arc/arc-*/*.html"))
    for file in html_files:
        report.analysis_files = [file]
        report.search_files(["cellranger_arc"])

    config.preserve_module_raw_data = True
    m = MultiqcModule()
    assert m.saved_raw_data is not None
    assert len(m.saved_raw_data) > 0

    data = m.saved_raw_data["multiqc_cellranger_arc"]

    # Check at least one sample exists
    assert len(data) > 0, "No samples found in the data"

    # Test each sample for proper ATAC/GEX separation
    for sample_name, sample_data in data.items():
        # Verify ATAC metrics exist with namespace prefix
        assert "ATAC_Sequenced read pairs" in sample_data, (
            f"Sample {sample_name}: Missing 'ATAC_Sequenced read pairs' (namespace not applied)"
        )
        assert "ATAC_Valid barcodes" in sample_data, (
            f"Sample {sample_name}: Missing 'ATAC_Valid barcodes' (namespace not applied)"
        )

        # Verify GEX metrics exist with namespace prefix
        assert "GEX_Sequenced read pairs" in sample_data, (
            f"Sample {sample_name}: Missing 'GEX_Sequenced read pairs' (namespace not applied)"
        )
        assert "GEX_Valid barcodes" in sample_data, (
            f"Sample {sample_name}: Missing 'GEX_Valid barcodes' (namespace not applied)"
        )

        # Critical test: Verify ATAC and GEX values are different (the bug was they collided)
        atac_reads = sample_data["ATAC_Sequenced read pairs"]
        gex_reads = sample_data["GEX_Sequenced read pairs"]
        assert atac_reads != gex_reads, (
            f"Sample {sample_name}: ATAC and GEX read pairs should be different! "
            f"ATAC={atac_reads}, GEX={gex_reads}. This indicates namespace collision."
        )

        atac_barcodes = sample_data["ATAC_Valid barcodes"]
        gex_barcodes = sample_data["GEX_Valid barcodes"]
        assert atac_barcodes != gex_barcodes, (
            f"Sample {sample_name}: ATAC and GEX valid barcodes should be different! "
            f"ATAC={atac_barcodes}, GEX={gex_barcodes}. This indicates namespace collision."
        )

        # Verify values are numeric
        assert isinstance(atac_reads, (int, float)), f"Sample {sample_name}: ATAC read pairs should be numeric"
        assert isinstance(gex_reads, (int, float)), f"Sample {sample_name}: GEX read pairs should be numeric"
        assert isinstance(atac_barcodes, (int, float)), f"Sample {sample_name}: ATAC valid barcodes should be numeric"
        assert isinstance(gex_barcodes, (int, float)), f"Sample {sample_name}: GEX valid barcodes should be numeric"

        # Log the values for verification
        print(f"\n{sample_name} metrics:")
        print(f"  ATAC Sequenced read pairs: {atac_reads:,.0f}")
        print(f"  GEX Sequenced read pairs: {gex_reads:,.0f}")
        print(f"  ATAC Valid barcodes: {atac_barcodes}")
        print(f"  GEX Valid barcodes: {gex_barcodes}")
