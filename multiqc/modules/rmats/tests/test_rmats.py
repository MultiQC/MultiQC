import pytest

from multiqc import config, report, reset
from multiqc.modules.rmats import MultiqcModule
from multiqc.modules.rmats.read_outcomes import parse_read_outcomes_file
from multiqc.modules.rmats.summary import parse_summary_file, summary_name_and_cutoff

OUTCOME_KEYS = [
    "USED",
    "NOT_PAIRED",
    "NOT_NH_1",
    "NOT_EXPECTED_CIGAR",
    "NOT_EXPECTED_READ_LENGTH",
    "NOT_EXPECTED_STRAND",
    "EXON_NOT_MATCHED_TO_ANNOTATION",
    "JUNCTION_NOT_MATCHED_TO_ANNOTATION",
    "CLIPPED",
]


def outcome_block(bam: str, used: int = 90, read_length: int = 10, clipped: int = 0) -> str:
    counts = dict.fromkeys(OUTCOME_KEYS, 0)
    counts["USED"] = used
    counts["NOT_EXPECTED_READ_LENGTH"] = read_length
    counts["CLIPPED"] = clipped
    lines = [bam] + [f"{k}: {v}" for k, v in counts.items()] + [f"TOTAL_FOR_BAM: {sum(counts.values())}"]
    return "\n".join(lines) + "\n"


SUMMARY_V42 = """EventType\tEventTypeDescription\tTotalEventsJC\tTotalEventsJCEC\tSignificantEventsJC\tSigEventsJCSample1HigherInclusion\tSigEventsJCSample2HigherInclusion\tSignificantEventsJCEC\tSigEventsJCECSample1HigherInclusion\tSigEventsJCECSample2HigherInclusion
SE\tskipped exon\t21\t23\t1\t1\t0\t1\t1\t0
A5SS\talternative 5' splice sites\t7\t8\t0\t0\t0\t0\t0\t0
A3SS\talternative 3' splice sites\t1\t1\t1\t0\t1\t1\t0\t1
MXE\tmutually exclusive exons\t1\t3\t0\t0\t0\t0\t0\t0
RI\tretained intron\t9\t11\t0\t0\t0\t0\t0\t0
"""

SUMMARY_V41 = """EventType\tTotalEventsJC\tTotalEventsJCEC\tSignificantEventsJC\tSigEventsJCSample1HigherInclusion\tSigEventsJCSample2HigherInclusion\tSignificantEventsJCEC\tSigEventsJCECSample1HigherInclusion\tSigEventsJCECSample2HigherInclusion
SE\t21\t23\t1\t1\t0\t1\t1\t0
A5SS\t7\t8\t0\t0\t0\t0\t0\t0
A3SS\t1\t1\t1\t0\t1\t1\t0\t1
MXE\t1\t3\t0\t0\t0\t0\t0\t0
RI\t9\t11\t0\t0\t0\t0\t0\t0
"""


class TestParseReadOutcomes:
    def test_multiple_bams(self):
        text = outcome_block("/data/a.bam") + outcome_block("/data/b.bam", used=50, read_length=50)
        result = parse_read_outcomes_file(text.splitlines(True))
        assert list(result) == ["/data/a.bam", "/data/b.bam"]
        assert result["/data/b.bam"]["USED"] == 50
        assert result["/data/b.bam"]["TOTAL_FOR_BAM"] == 100

    def test_windows_line_endings(self):
        text = outcome_block("C:\\data\\a.bam").replace("\n", "\r\n")
        result = parse_read_outcomes_file(text.splitlines(True))
        assert result["C:\\data\\a.bam"]["USED"] == 90

    def test_missing_outcome_raises(self):
        text = outcome_block("a.bam").replace("CLIPPED: 0\n", "")
        with pytest.raises(ValueError, match="missing read outcomes: CLIPPED"):
            parse_read_outcomes_file(text.splitlines(True))

    def test_total_mismatch_raises(self):
        text = outcome_block("a.bam").replace("TOTAL_FOR_BAM: 100", "TOTAL_FOR_BAM: 101")
        with pytest.raises(ValueError, match="sum to 100"):
            parse_read_outcomes_file(text.splitlines(True))

    def test_truncated_raises(self):
        text = "".join(outcome_block("a.bam").splitlines(True)[:-1])
        with pytest.raises(ValueError, match="no TOTAL_FOR_BAM"):
            parse_read_outcomes_file(text.splitlines(True))

    def test_unknown_outcome_is_kept(self):
        text = outcome_block("a.bam").replace("CLIPPED: 0\n", "CLIPPED: 0\nNEW_REASON: 5\n")
        text = text.replace("TOTAL_FOR_BAM: 100", "TOTAL_FOR_BAM: 105")
        result = parse_read_outcomes_file(text.splitlines(True))
        assert result["a.bam"]["NEW_REASON"] == 5


class TestParseSummary:
    @pytest.mark.parametrize("text", [SUMMARY_V42, SUMMARY_V41])
    def test_versions(self, text):
        result = parse_summary_file(text.splitlines(True))
        assert result["SE_TotalEventsJCEC"] == 23
        assert result["A3SS_SigEventsJCSample2HigherInclusion"] == 1
        assert len(result) == 5 * 8

    def test_missing_column_raises(self):
        text = SUMMARY_V41.replace("\tSigEventsJCECSample2HigherInclusion", "")
        with pytest.raises(ValueError, match="missing columns"):
            parse_summary_file(text.splitlines(True))

    def test_missing_event_type_raises(self):
        text = "".join(line for line in SUMMARY_V41.splitlines(True) if not line.startswith("RI\t"))
        with pytest.raises(ValueError, match="RI"):
            parse_summary_file(text.splitlines(True))

    def test_non_integer_raises(self):
        text = SUMMARY_V41.replace("SE\t21", "SE\tNA")
        with pytest.raises(ValueError):
            parse_summary_file(text.splitlines(True))


class TestSummaryName:
    @pytest.mark.parametrize(
        "fn,root,expected",
        [
            ("summary.txt", "results/PC3E-GS689", ("PC3E-GS689", None)),
            ("summary.txt", "", ("rmats", None)),
            (
                "PC3E-GS689_FDR_0.01_IncLevelDifference_0.1.txt",
                "out",
                ("PC3E-GS689_FDR_0.01_IncLevelDifference_0.1", "FDR <= 0.01, |IncLevelDifference| >= 0.1"),
            ),
            (
                "summary_PValue_0.05_IncLevelDifference_0.0.txt",
                "od",
                ("od_PValue_0.05_IncLevelDifference_0.0", "PValue <= 0.05, |IncLevelDifference| >= 0.0"),
            ),
        ],
    )
    def test_names(self, fn, root, expected):
        assert summary_name_and_cutoff(fn, root) == expected


@pytest.fixture(autouse=True)
def reset_multiqc():
    reset()
    config.preserve_module_raw_data = True
    yield
    reset()


def run_module(tmp_path, files: dict[str, str]):
    for rel, text in files.items():
        path = tmp_path / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    report.analysis_files = [tmp_path]
    report.search_files(["rmats"])
    return MultiqcModule()


class TestModule:
    def test_sample_names_and_alerts(self, tmp_path):
        text = (
            outcome_block("/tmp/x/sampleA.bam", used=99, read_length=1)
            + outcome_block("/tmp/2024_01_01_bam1_1/Aligned.sortedByCoord.out.bam", used=20, read_length=80)
            + outcome_block("C:\\data\\sampleC.bam", used=0, read_length=0, clipped=5)
        )
        m = run_module(tmp_path, {"2024-01-01-00_00_00_read_outcomes_by_bam.txt": text})
        names = set(m.saved_raw_data["multiqc_rmats_read_outcomes"])
        assert names == {"sampleA", "2024_01_01_bam1_1", "sampleC"}
        section = m.sections[0]
        alerts = {(a.level, "--allow-clipping" in a.message): a.affected_samples for a in section.alerts}
        assert alerts == {
            ("danger", False): ["sampleC"],
            ("warning", False): ["2024_01_01_bam1_1"],
            ("warning", True): ["sampleC"],
        }

    def test_summary_only(self, tmp_path):
        m = run_module(
            tmp_path,
            {"PC3E-GS689/summary.txt": SUMMARY_V42, "out/GS689-PC3E.rmats_summary.txt": SUMMARY_V41},
        )
        assert set(m.saved_raw_data["multiqc_rmats_summary"]) == {"PC3E-GS689", "GS689-PC3E"}
        assert report.general_stats_data == {}
        assert [s.anchor for s in m.sections] == ["rmats_event_counts", "rmats_significant_events"]

    def test_no_significant_events_keeps_section(self, tmp_path):
        zero = "\n".join(
            "\t".join(cols[:4] + ["0"] * 6) if i else line
            for i, line in enumerate(SUMMARY_V42.splitlines())
            for cols in [line.split("\t")]
        )
        m = run_module(tmp_path, {"c1/summary.txt": zero + "\n"})
        significant = next(s for s in m.sections if s.anchor == "rmats_significant_events")
        assert significant.plot_anchor is None
        assert significant.alerts

    def test_skips_malformed_file(self, tmp_path):
        bad = SUMMARY_V41.replace("SE\t21", "SE\tNA")
        m = run_module(tmp_path, {"a/summary.txt": SUMMARY_V41, "b/summary.txt": bad})
        assert set(m.saved_raw_data["multiqc_rmats_summary"]) == {"a"}
