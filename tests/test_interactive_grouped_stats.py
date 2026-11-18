import multiqc


def test_get_general_stats_grouped_samples(data_dir, tmp_path):
    cfg = tmp_path / "multiqc_config.yaml"
    cfg.write_text('table_sample_merge:\n  R1: "_1"\n  R2: "_2"\n')
    multiqc.parse_logs(
        data_dir / "modules/fastqc/groups/test_1_fastqc.zip",
        data_dir / "modules/fastqc/groups/test_2_fastqc.zip",
        data_dir / "modules/fastqc/groups/single_fastqc.zip",
        config_files=[cfg],
        run_modules=["fastqc"],
        no_version_check=True,
        strict=True,
    )
    # Raw FastQC reports each contain 266736 reads, with 50% and 48% GC.
    stats = multiqc.get_general_stats_data()
    assert stats["test"]["FastQC.total_sequences"] == 533472
    assert stats["test"]["FastQC.percent_gc"] == 49
    assert stats["test R1"]["FastQC.total_sequences"] == 266736
    assert stats["test R1"]["FastQC.percent_gc"] == 50
    assert stats["test R2"]["FastQC.total_sequences"] == 266736
    assert stats["test R2"]["FastQC.percent_gc"] == 48
    assert stats["single"]["FastQC.total_sequences"] == 266736
    assert stats["single"]["FastQC.percent_gc"] == 50
    for sample, data in stats.items():
        assert multiqc.get_general_stats_data(sample=sample) == data
    assert multiqc.get_general_stats_data(sample="missing") == {}
