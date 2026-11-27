from multiqc import report

GROUPING_TSV = (
    "accepted_sam_records\tdiscarded_non_pf\tdiscarded_poor_alignment\tdiscarded_ns_in_umi\tdiscarded_umis_to_short\n"
    "900\t10\t50\t30\t10\n"
)
POSITION_TSV = (
    "position_group_size\tcount\tfraction\tfraction_gt_or_eq_position_group_size\n1\t80\t0.8\t1\n2\t20\t0.2\t0.2\n"
)


def test_grouping_metrics_bar_and_position_sizes(run_fgumi):
    module = run_fgumi({"S1.grouping_metrics.txt": GROUPING_TSV, "S1.position_group_sizes.txt": POSITION_TSV})
    assert module.samples_parsed_by_tool["grouping"] == {"S1"}
    assert module.saved_raw_data["multiqc_fgumi_grouping_metrics"]["S1"] == {
        "accepted_sam_records": 900,
        "discarded_non_pf": 10,
        "discarded_poor_alignment": 50,
        "discarded_ns_in_umi": 30,
        "discarded_umis_to_short": 10,
    }
    assert module.saved_raw_data["multiqc_fgumi_position_group_sizes"]["S1"] == {1: 80, 2: 20}
    assert {"fgumi_grouping_metrics", "fgumi_position_group_sizes"} <= set(report.plot_by_id)
