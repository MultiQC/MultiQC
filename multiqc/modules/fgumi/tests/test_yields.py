from .conftest import general_stats, line_points

DUPLEX_YIELD = (
    "fraction\tread_pairs\tcs_families\tss_families\tds_families\tds_duplexes\tds_fraction_duplexes"
    "\tds_fraction_duplexes_ideal\n"
    "0.5\t500\t200\t150\t100\t40\t0.4\t0.5\n"
    "1\t1000\t350\t300\t200\t100\t0.5\t0\n"
)
SIMPLEX_YIELD = (
    "fraction\tread_pairs\tcs_families\tss_families\tmean_ss_family_size\tss_singletons\tss_singleton_fraction"
    "\tss_consensus_families\n"
    "1\t1000\t300\t250\t4\t25\t0.1\t200\n"
)


def test_duplex_yield_curves_and_ratio_guards_zero_ideal(run_fgumi):
    module = run_fgumi({"B.duplex_yield_metrics.txt": DUPLEX_YIELD})
    raw = module.saved_raw_data["multiqc_fgumi_duplex_yield"]["B"]
    assert (raw["duplexes_500"], raw["duplexes_1000"]) == (40, 100)
    assert line_points("fgumi_duplex_yield", 1, "B") == {500: 50.0, 1000: 0.0}
    # actual/ideal: 0.4/0.5 = 0.8 at 500 pairs; ideal 0 at 1000 pairs gives no point rather than inf
    assert line_points("fgumi_duplex_yield", 3, "B") == {500: 0.8}
    assert "ratio_1000" not in raw
    # General Statistics use the full-depth row, not the first (downsampled) one.
    assert general_stats("B")["ds_duplexes"] == 100


def test_simplex_yield(run_fgumi):
    module = run_fgumi({"A.simplex_yield_metrics.txt": SIMPLEX_YIELD})
    assert module.saved_raw_data["multiqc_fgumi_simplex_yield"]["A"]["consensus_families_1000"] == 200
    assert general_stats("A")["ss_consensus_families"] == 200


def test_non_finite_ideal_fraction_gives_no_ratio_point_not_zero(run_fgumi):
    rows = DUPLEX_YIELD.replace("0.5\t500\t200\t150\t100\t40\t0.4\t0.5", "0.5\t500\t200\t150\t100\t40\t0.4\tInfinity")
    module = run_fgumi({"B.duplex_yield_metrics.txt": rows})
    # Neither depth has a finite ratio, so the ratio is left out of the data file (and the plot has no points).
    assert not [
        column for column in module.saved_raw_data["multiqc_fgumi_duplex_yield"]["B"] if column.startswith("ratio")
    ]
