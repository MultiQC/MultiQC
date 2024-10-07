from multiqc import report
from multiqc.modules.umicollapse import MultiqcModule
from multiqc.plots.table_object import InputRow


def test_parse(tmp_path):
    # File without a file name or content match.
    f1 = tmp_path / "SAMPLE.log"
    f1.write_text("Irrelevant file")

    # File without file name match.
    f2 = tmp_path / "SAMPLE.log"
    f2.write_text("UMI collapsing finished in 1077.717 seconds!")

    # File with match.
    f3 = tmp_path / "SRR19887568.umi_dedup.sorted_UMICollapse.log"
    f3.write_text("""\
Arguments	[bam, -i, SRR19887568.sorted.bam, -o, SRR19887568.umi_dedup.sorted.bam, --paired, --two-pass, --remove-unpaired, --remove-chimeric, --algo, dir]
Number of input reads	53490614
Number of removed unmapped reads	24260
Number of unpaired reads	0
Number of chimeric reads	0
Number of unremoved reads	53466354
Number of unique alignment positions	21532282
Average number of UMIs per alignment position	1.3478738110526325
Max number of UMIs over all alignment positions	165
Number of reads after deduplicating	27208411
UMI collapsing finished in 1077.717 seconds!
""")

    report.analysis_files = [f1, f2, f3]
    report.search_files(["umicollapse"])
    m = MultiqcModule()
    assert len(m.saved_raw_data) == 1
    print(m.saved_raw_data)

    data = m.saved_raw_data["multiqc_umicollapse"]

    keys = data.keys()
    assert len(keys) == 1
    assert list(keys)[0] == "SRR19887568"

    assert len(m.sections) == 2
    assert m.sections[0].name == "Deduplicated Reads"
    assert m.sections[1].name == "UMI Stats"

    assert len(report.general_stats_data) > 0
    assert report.general_stats_data[-1] == {
        "SRR19887568": [
            InputRow(
                sample="SRR19887568",
                data={
                    "input_reads": 53490614,
                    "dedup_input_reads": 53466354,
                    "positions_deduplicated": 21532282,
                    "mean_umi_per_pos": 1.3478738110526325,
                    "max_umi_per_pos": 165,
                    "dedup_output_reads": 27208411,
                    "dedup_percent_passing": 50.89,
                    "removed_reads": 24260,
                    "dedup_removed_reads": 26257943,
                },
            )
        ]
    }
