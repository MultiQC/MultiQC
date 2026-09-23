import pytest

CLIP = (
    "read_type\treads\treads_unmapped\treads_clipped_pre\treads_clipped_post\treads_clipped_five_prime"
    "\treads_clipped_three_prime\treads_clipped_overlapping\treads_clipped_extending\tbases\tbases_clipped_pre"
    "\tbases_clipped_post\tbases_clipped_five_prime\tbases_clipped_three_prime\tbases_clipped_overlapping"
    "\tbases_clipped_extending\n"
    "Fragment\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\t0\n"
    "All\t100\t0\t0\t30\t5\t10\t12\t3\t10000\t0\t400\t50\t100\t200\t50\n"
)
FILTER = "total_reads\tpassed_reads\tfailed_reads\tpass_rate\n1000\t900\t100\t0.9\n"
COPY_UMI = "total_records\trx_written\trx_overwritten\tnames_trimmed\n100\t100\t4\t100\n"
RETAG = "operation\tkind\trecords_applied\tdst_overwritten\tsrc_missing\nRX::copy::BX\tcopy\t90\t0\t10\n"
DOWNSAMPLE = "family_size\tcount\n1\t40\n2\t10\n"
FAMILY_SIZES = "family_size\tcount\tfraction\tfraction_gt_or_eq_family_size\n1\t60\t0.6\t1\n"


def test_small_command_metrics(run_fgumi):
    # One file per tool; names end in plain ".txt" so MultiQC's default cleaning yields exactly these samples.
    module = run_fgumi(
        {"Clip1.txt": CLIP, "Filt1.txt": FILTER, "Copy1.txt": COPY_UMI, "Retag1.txt": RETAG, "Down1.txt": DOWNSAMPLE}
    )
    raw = module.saved_raw_data
    assert raw["multiqc_fgumi_clip"]["Clip1"]["reads"] == {
        "five_prime": 5,
        "three_prime": 10,
        "overlapping": 12,
        "extending": 3,
    }
    filt = raw["multiqc_fgumi_filter"]["Filt1"]
    assert (filt["passed"], filt["failed"]) == (900, 100)
    assert filt["pass_rate"] == pytest.approx(90.0)
    assert raw["multiqc_fgumi_copy_umi"]["Copy1"] == {"rx_written": 100, "rx_overwritten": 4, "names_trimmed": 100}
    assert raw["multiqc_fgumi_retag"]["Retag1 (RX::copy::BX)"]["records_applied"] == 90
    assert raw["multiqc_fgumi_downsample"]["Down1"] == {1: 40, 2: 10}


def test_downsample_pattern_does_not_steal_family_size_histograms(run_fgumi):
    module = run_fgumi({"S1.family_sizes.txt": FAMILY_SIZES})
    assert module.samples_parsed_by_tool["family_sizes"] == {"S1"}
    assert module.samples_parsed_by_tool["commands"] == set()


LEGACY_FILTER = "total_reads\t1000\npassed_reads\t900\nfailed_reads\t100\npass_rate\t0.9000\n"


def test_filter_stats_legacy_headerless_layout(run_fgumi):
    # fgumi <= 0.7.0 writes filter --stats as headerless key/value rows.
    module = run_fgumi({"Filt1.txt": LEGACY_FILTER})
    filt = module.saved_raw_data["multiqc_fgumi_filter"]["Filt1"]
    assert (filt["passed"], filt["failed"]) == (900, 100)
    assert filt["pass_rate"] == pytest.approx(90.0)


def test_clip_file_without_all_row_warns(run_fgumi, caplog):
    no_all = "\n".join(line for line in CLIP.splitlines() if not line.startswith("All")) + "\n"
    run_fgumi({"Clip1.txt": no_all, "Filt1.txt": FILTER})
    assert any("Clip1.txt" in r.message and "All" in r.message for r in caplog.records)


def test_retag_honours_sample_filters_on_the_sample_name(run_fgumi):
    from multiqc import config

    original = config.sample_names_only_include
    config.sample_names_only_include = ["Retag1"]
    try:
        module = run_fgumi({"Retag1.txt": RETAG, "Filt1.txt": FILTER})
    finally:
        config.sample_names_only_include = original
    assert module.samples_parsed_by_tool["commands"] == {"Retag1"}
    assert list(module.saved_raw_data["multiqc_fgumi_retag"]) == ["Retag1 (RX::copy::BX)"]
