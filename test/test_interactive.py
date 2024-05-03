import multiqc


def test_parse_logs():
    multiqc.parse_logs(
        "../test_data/data/modules/fastp/SAMPLE.json",
        "../test_data/data/modules/fastp/single_end",
        extra_fn_clean_exts=["_1", "_S10_R1_001"],
    )
    assert multiqc.list_samples() == ["SRR5442949", "smalltest"]
    assert multiqc.list_modules() == ["fastp"]

    multiqc.parse_logs(
        "../test_data/data/modules/quast/full_metaquast_run",
        ignore_samples=["meta_contigs_2"],
    )

    assert multiqc.list_samples() == ["SRR5442949", "meta_contigs_1", "smalltest"]
    assert multiqc.list_modules() == ["fastp", "QUAST"]
