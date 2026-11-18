import multiqc
from multiqc import report


def test_parquet_merge_preserves_section_order(data_dir, tmp_path):
    sources = [data_dir / "modules/fastp/SAMPLE.json", data_dir / "modules/fastp/single_end"]
    parquets = []
    for i, source in enumerate(sources):
        multiqc.parse_logs(source, run_modules=["fastp"], no_version_check=True, strict=True)
        if i == 0:
            expected_sections = [s.anchor for s in report.modules[0].sections]
        out = tmp_path / str(i)
        multiqc.write_report(output_dir=out, make_report=False, no_version_check=True)
        parquets.append(out / "multiqc_data/multiqc.parquet")
        multiqc.reset()

    # Single-file reload is the control; merging must keep the same section sequence.
    multiqc.parse_logs(parquets[0], no_version_check=True, strict=True)
    assert [s.anchor for s in report.modules[0].sections] == expected_sections
    multiqc.reset()
    multiqc.parse_logs(*parquets, no_version_check=True, strict=True)
    merged_sections = [s.anchor for s in report.modules[0].sections]
    assert merged_sections == expected_sections
    multiqc.reset()
    multiqc.parse_logs(*sources, run_modules=["fastp"], no_version_check=True, strict=True)
    assert [s.anchor for s in report.modules[0].sections] == merged_sections
