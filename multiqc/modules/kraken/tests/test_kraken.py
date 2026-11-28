from io import StringIO
from pathlib import Path

from multiqc import parse_logs, report
from multiqc.modules.kraken.kraken import parse_logs as parse_kraken_report
from multiqc.types import Anchor

KRAKEN_REPORT = (
    " 25.00\t1\t1\t0\t0\tU\t0\tunclassified\n"
    " 75.00\t3\t0\t57\t8\tR\t1\troot\n"
    " 75.00\t3\t0\t57\t8\tG\t10\t  Test genus\n"
    " 75.00\t3\t3\t57\t8\tS\t11\t    Test genus observed\n"
)
ZERO_COUNT_SPECIES = "  0.00\t0\t0\t0\t0\tS\t12\t    Test genus absent\n"


def test_zero_distinct_minimizers_preserve_species_counts():
    counts, duplication = parse_kraken_report({"fn": "padded.txt", "f": StringIO(KRAKEN_REPORT + ZERO_COUNT_SPECIES)})

    assert counts["S"] == {"Test genus observed": 3, "Test genus absent": 0}
    assert duplication == {"Test genus observed": 57 / 8}


def test_zero_count_kraken_sample_does_not_drop_other_samples(tmp_path: Path):
    (tmp_path / "positive.txt").write_text(KRAKEN_REPORT)
    (tmp_path / "padded.txt").write_text(KRAKEN_REPORT + ZERO_COUNT_SPECIES)

    parse_logs(tmp_path, run_modules=["kraken"], strict=True, no_version_check=True)

    counts = report.saved_raw_data["multiqc_kraken"]
    assert set(counts) == {"positive", "padded"}
    assert counts["positive"]["S"] == {"Test genus observed": 3}
    assert counts["padded"]["S"] == {"Test genus observed": 3, "Test genus absent": 0}
    assert Anchor("kraken-top-duplication_plot") in report.plot_by_id
