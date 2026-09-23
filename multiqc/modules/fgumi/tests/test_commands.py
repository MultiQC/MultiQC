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
